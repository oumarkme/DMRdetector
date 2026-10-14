"""Smoothing with the FFT of each chromosome in a child process.

The bytes equal those of in-process FFTs, in the fast path and in the pool path.
"""
import hashlib
import os
import sys

import numpy as np
import pytest


sys.path.insert(0, os.path.dirname(__file__))
from dmrdetector.fileio import bed as dio  # noqa: E402
from dmrdetector.steps import atlas  # noqa: E402
from dmrdetector.util import fftsmooth  # noqa: E402


def _chrom(seed, span=60000, n=700):
    rng = np.random.default_rng(seed)
    start = np.unique(rng.integers(1000, span, n))
    return (start, start + 1, rng.integers(5, 60, start.size).astype(float),
            rng.random(start.size), rng.random(start.size) * 0.2)


@pytest.mark.parametrize("inner", [1, 3])
def test_child_equals_in_process(inner):
    start, end, cov, raw, phi = _chrom(1)
    k = atlas.gaussian_kernel(50)
    a = atlas.smooth_chromosome(start, end, cov, raw, phi, k, inner)
    b = atlas.smooth_chromosome(start, end, cov, raw, phi, k, inner, child=True)
    assert (a[0] == b[0]).all() and (a[1] == b[1]).all()  # bitwise
    assert a[0].dtype == b[0].dtype == np.float64 and a[0].shape == b[0].shape


def test_child_failure_is_an_error(monkeypatch):
    start, end, cov, raw, phi = _chrom(2)
    monkeypatch.setattr(sys, "executable", "/bin/false")
    with pytest.raises(RuntimeError, match="FFT child"):
        atlas.smooth_chromosome(start,
                                end,
                                cov,
                                raw,
                                phi,
                                atlas.gaussian_kernel(50),
                                1,
                                child=True)


def _raw_atlas(path):
    rows = []
    for ci, (c, span) in enumerate(
        (("chr1", 90000), ("chr2", 70000), ("chr3", 50000), ("chrX", 40000),
         ("chrY", 30000))):
        start, _, cov, raw, phi = _chrom(10 + ci, span, 500)
        rows += [(c, s, s + 1, 5, int(v), float(p), float(m))
                 for s, v, p, m in zip(start, cov, phi, raw)]
    dio.write_bgzip_text(path, [
        dio.join_columns(
            [[r[i] if i == 0 else str(r[i]) for r in rows] for i in range(7)])
    ])


def _md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("threads", [2, 9])
def test_smooth_atlas_same_bytes_with_and_without_child(tmp_path, monkeypatch,
                                                        legacy, threads):
    raw = str(tmp_path / "raw.bed.gz")
    _raw_atlas(raw)
    if legacy:
        monkeypatch.setattr(atlas, "_load_raw_atlas", lambda *a, **k: None)
    outs = {}
    for child in (False, True):
        out = str(tmp_path / f"ref_{child}.bed.gz")
        atlas.smooth_atlas(raw,
                           40,
                           out,
                           True,
                           threads=threads,
                           log=lambda m: None,
                           fft_child=child)
        outs[child] = out
    assert _md5(outs[True]) == _md5(outs[False])
    assert _md5(outs[True] + ".tbi") == _md5(outs[False] + ".tbi")


def test_module_runs_as_child_with_log(tmp_path, monkeypatch, capfd):
    monkeypatch.setenv("DMRDETECTOR_FFT_LOG", "1")
    start, end, cov, raw, phi = _chrom(3)
    atlas.smooth_chromosome(start,
                            end,
                            cov,
                            raw,
                            phi,
                            atlas.gaussian_kernel(50),
                            1,
                            child=True)
    assert "fft child: length" in capfd.readouterr().err
    assert fftsmooth.__name__ == "dmrdetector.util.fftsmooth"


def test_memory_estimate_fits_the_minimum_budget(tmp_path):
    from dmrdetector import build
    from dmrdetector import runtime
    chr1 = runtime.SCOPE_LARGEST_SPAN["autosomal"]
    assert 9.9 <= runtime.smooth_gb(
        chr1
    ) <= 10.6  # child peak measured at 9.88 GB (250 million padded bases)
    assert 17.7 <= runtime.smooth_gb(chr1,
                                     3) <= 18.5  # 17.66 GB with 3 FFT threads
    assert runtime.smooth_gb(chr1 // 2) < runtime.smooth_gb(chr1)
    p = tmp_path / "s_autosomal.bed.gz"
    with open(p, "wb") as fh:
        fh.truncate(int(2.51e9))  # size of a real raw atlas
    assert build.scope_need_gb(
        [str(p)], "autosomal") <= runtime.MIN_THREADS * runtime.GB_PER_THREAD


def test_child_ignores_the_working_directory(tmp_path, monkeypatch):
    """Checks that the child does not import from the working directory.

    A numpy.py or scipy/ there must not be imported (python -P).
    """
    marker = tmp_path / "decoy_imported"
    (tmp_path / "numpy.py").write_text(
        f"open({str(marker)!r}, 'w').close()\nraise ImportError('decoy "
        "numpy')\n")
    (tmp_path / "scipy").mkdir()
    (tmp_path / "scipy" / "__init__.py").write_text(
        f"open({str(marker)!r}, 'w').close()\nraise ImportError('decoy "
        "scipy')\n")
    start, end, cov, raw, phi = _chrom(3)
    k = atlas.gaussian_kernel(50)
    a = atlas.smooth_chromosome(start, end, cov, raw, phi, k, 1)
    monkeypatch.chdir(tmp_path)
    b = atlas.smooth_chromosome(start, end, cov, raw, phi, k, 1, child=True)
    assert (a[0] == b[0]).all() and (a[1] == b[1]).all()
    assert not marker.exists()


def test_multiply_complex_is_the_plain_formula_in_any_block_size():
    rng = np.random.default_rng(5)
    a = rng.normal(size=4001) + 1j * rng.normal(size=4001)
    b = rng.normal(size=4001) + 1j * rng.normal(size=4001)
    want = np.array([
        complex(x.real * y.real - x.imag * y.imag,
                x.real * y.imag + x.imag * y.real) for x, y in zip(a, b)
    ])  # Python floats: every product and sum is rounded on its own
    for chunk in (1, 7, 1 << 20):
        c = a.copy()
        fftsmooth.multiply_complex(c, b, chunk)
        assert (c.view(np.float64) == want.view(np.float64)).all(), chunk
    np.testing.assert_allclose(a * b, want, rtol=1e-15)
    with pytest.raises(ValueError):
        fftsmooth.multiply_complex(a.copy(), b[:-1])
    with pytest.raises(ValueError):
        fftsmooth.multiply_complex(a.copy(), b.astype(np.complex64))
