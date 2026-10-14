"""FFT convolution of one chromosome, in this process or in a child process.

scipy.fft keeps the plan of every transform length it has seen (about 8 bytes
per padded base) until the process ends, and chromosomes differ in length, so a
process that smooths several holds all their plans. A child that smooths one
chromosome and exits returns that memory. The child runs this module (`python -m
dmrdetector.util.fftsmooth`) and gets its input over a pipe; the function and
its arguments are the same, so the numbers are identical. No fork, start method
or shared memory is involved. DMRDETECTOR_FFT_LOG=1 makes the child report its
peak memory on stderr.
"""

from __future__ import annotations

import os
import pickle
import subprocess
import sys

import numpy as np


def multiply_complex(a: np.ndarray,
                     b: np.ndarray,
                     chunk: int = 1 << 20) -> None:
    """Sets a = a * b for complex128 arrays of one shape, without fused ops.

    The result is (ar * br - ai * bi) + i (ar * bi + ai * br), with every
    product and sum rounded to float64 on its own. `a *= b` of numpy does the
    same arithmetic in one loop, and numpy 2.5.3 runs that loop with fused
    multiply-add instructions on a CPU with AVX2 and FMA3 and with separate
    multiplications and additions on another (X86_V3 dispatch), which changes
    the last bit of some products. Separate ufunc calls cannot be fused, so the
    result is the same on every CPU. The work goes in blocks so that the
    temporary arrays stay small.

    Args:
      a: Complex128 array, changed in place.
      b: Complex128 array of the same shape.
      chunk: Number of elements per block.
    """
    if (a.shape != b.shape or a.dtype != np.complex128 or
            b.dtype != np.complex128 or not a.flags.c_contiguous or
            not b.flags.c_contiguous):
        raise ValueError(
            "multiply_complex needs contiguous complex128 arrays of one shape")
    af = a.view(np.float64).reshape(-1, 2)
    bf = b.view(np.float64).reshape(-1, 2)
    for lo in range(0, af.shape[0], chunk):
        ar, ai = af[lo:lo + chunk, 0], af[lo:lo + chunk, 1]
        br, bi = bf[lo:lo + chunk, 0], bf[lo:lo + chunk, 1]
        re = ar * br
        re -= ai * bi
        im = ar * bi
        im += ai * br
        af[lo:lo + chunk, 0] = re
        af[lo:lo + chunk, 1] = im


def smooth_core(idx: np.ndarray,
                num: np.ndarray,
                den: np.ndarray,
                kernel: np.ndarray,
                fft_len: int,
                span: int,
                inner: int = 1) -> np.ndarray:
    """Returns REF_METH = conv(kernel, num) / conv(kernel, den) at the sites.

    The sites are `idx`. One zero-padded buffer is reused for the kernel and
    both signals and the product is formed in place (about 3 complex-half
    arrays, 24 bytes per padded base). With inner > 1 the independent transforms
    (kernel, the two signals, then the two inverses) run in up to 3 threads
    (about 40 bytes per padded base at the peak); each transform is the same
    single-threaded call, so the result is bit-identical to inner = 1.
    """
    from scipy import fft as sfft
    k_len = kernel.size
    radius = k_len // 2

    def fwd_kernel():
        """Returns the real FFT of the zero-padded kernel."""
        buf = np.zeros(fft_len)
        buf[:k_len] = kernel
        return sfft.rfft(buf, n=fft_len, overwrite_x=True)

    def fwd(vals):
        """Returns the real FFT of `vals` placed at the sites `idx`."""
        buf = np.zeros(fft_len)
        buf[idx] = vals
        return sfft.rfft(buf, n=fft_len, overwrite_x=True)

    def inv(fv, kf):
        """Returns the inverse FFT of `fv * kf` at the sites `idx`."""
        multiply_complex(fv, kf)
        r = sfft.irfft(fv, n=fft_len, overwrite_x=True)
        del fv
        return r[radius:radius + span][idx]

    if inner > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(3) as ex:
            fk, f1, f2 = ex.submit(fwd_kernel), ex.submit(fwd, num), ex.submit(
                fwd, den)
            kf = fk.result()
            g1, g2 = ex.submit(inv, f1.result(),
                               kf), ex.submit(inv, f2.result(), kf)
            sm, sc = g1.result(), g2.result()
    else:
        kf = fwd_kernel()
        sm = inv(fwd(num), kf)
        sc = inv(fwd(den), kf)
    del kf
    ref = np.zeros(idx.size)
    ok = sc > 1e-9
    ref[ok] = sm[ok] / sc[ok]
    return ref


def _send(fh, obj) -> None:
    pickle.dump(obj, fh, protocol=4)


def _send_array(fh, a) -> None:
    fh.write(memoryview(np.ascontiguousarray(a)).cast("B"))


def _recv_array(fh, dtype, n):
    a = np.empty(n, dtype=dtype)
    view, got = memoryview(a).cast("B"), 0
    while got < view.nbytes:
        r = fh.readinto(view[got:])
        if not r:
            raise EOFError("pipe closed before the array was complete")
        got += r
    return a


def smooth_in_child(idx: np.ndarray,
                    num: np.ndarray,
                    den: np.ndarray,
                    kernel: np.ndarray,
                    fft_len: int,
                    span: int,
                    inner: int = 1) -> np.ndarray:
    """Runs `smooth_core` in a new Python process that exits afterwards.

    The result is the same array.
    """
    import dmrdetector
    root = os.path.dirname(
        os.path.dirname(os.path.abspath(dmrdetector.__file__)))
    env = dict(os.environ)
    env["PYTHONPATH"] = root + (os.pathsep + env["PYTHONPATH"]
                                if env.get("PYTHONPATH") else "")
    arrays = [np.ascontiguousarray(a) for a in (idx, num, den, kernel)]
    # -P: do not put the working directory first on sys.path (a stray numpy.py
    # or scipy/ there would be imported); the package location comes from
    # PYTHONPATH instead
    p = subprocess.Popen([sys.executable, "-P", "-m", __name__],
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         env=env)
    try:
        _send(
            p.stdin,
            dict(fft_len=int(fft_len),
                 span=int(span),
                 inner=int(inner),
                 arrays=[(a.dtype.str, a.size) for a in arrays]))
        for a in arrays:
            _send_array(p.stdin, a)
        p.stdin.close()
        head = pickle.load(p.stdout)
        ref = _recv_array(p.stdout, np.float64, head["n"])
    except (EOFError, BrokenPipeError, pickle.UnpicklingError) as e:
        p.kill()
        p.wait()
        raise RuntimeError(
            f"the FFT child process failed (exit status {p.returncode}): {e}"
        ) from e
    finally:
        if p.stdin and not p.stdin.closed:
            p.stdin.close()
        p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError(
            f"the FFT child process exited with status {p.returncode}")
    return ref


def _hwm_gb() -> float:
    try:
        with open("/proc/self/status") as fh:
            for ln in fh:
                if ln.startswith("VmHWM"):
                    return int(ln.split()[1]) / 1e6
    except OSError:
        pass
    return float("nan")


def _child_main() -> int:
    fin, fout = sys.stdin.buffer, sys.stdout.buffer
    head = pickle.load(fin)
    idx, num, den, kernel = (
        _recv_array(fin, np.dtype(d), n) for d, n in head["arrays"])
    ref = smooth_core(idx, num, den, kernel, head["fft_len"], head["span"],
                      head["inner"])
    if os.environ.get("DMRDETECTOR_FFT_LOG"):
        print(
            f"fft child: length {head['fft_len']}, inner {head['inner']}, "
            f"peak {_hwm_gb():.2f} GB",
            file=sys.stderr)
    _send(fout, dict(n=ref.size))
    _send_array(fout, ref)
    fout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(_child_main())
