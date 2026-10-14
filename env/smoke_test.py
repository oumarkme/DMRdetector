"""Image smoke test.

Checks the pinned interpreter, the locked versions, the imports, a numba jit
call, the compiled kernels of the image and a libdeflate round trip.

The kernels: a temporary NUMBA_CACHE_DIR is seeded from the cache baked into
the image (/opt/dmrdetector/numba-cache), a small data set runs through the
library (preprocess, raw atlas, smoothed atlas, CpG islands, so that every
module with baked entries calls a kernel), and the test asserts that kernels
were loaded from the cache and that none was compiled. Without the baked cache
the test prints "absent" and goes on, unless DMRDETECTOR_REQUIRE_CACHE=1 (set
by release builds: `--build-arg REQUIRE_CACHE=1` for Apptainer,
`--build-arg DMRDETECTOR_REQUIRE_CACHE=1` for Docker), when it fails.
"""
import atexit
import ctypes
import importlib.metadata as md
import importlib.util
import os
import platform
import re
import shutil
import sys
import tempfile
import zlib


LOCK = "/opt/dmrdetector/requirements.lock"
SEED_DIR = "/opt/dmrdetector/numba-cache"
REQUIRE_CACHE = os.environ.get("DMRDETECTOR_REQUIRE_CACHE") == "1"
WORK = tempfile.mkdtemp(prefix="dmrdetector-smoke-")
atexit.register(shutil.rmtree, WORK, ignore_errors=True)
# the folder the kernels are loaded from; numba reads the variable when it is
# imported. The program copies the baked cache into it at the first jit call.
os.environ["NUMBA_CACHE_DIR"] = os.path.join(WORK, "numba")
# the kernels of the worker processes would not show in this process
os.environ["DMRDETECTOR_EXECUTOR"] = "thread"

want_py = os.environ.get("PYTHON_VERSION")
assert want_py is None or platform.python_version() == want_py, (
    platform.python_version(), want_py)

# every locked distribution is installed at exactly the locked version
for line in open(LOCK):
    m = re.match(r"^([A-Za-z0-9._-]+)==(\S+)", line)
    if m:
        got = md.version(m.group(1))
        assert got == m.group(2), (m.group(1), got, m.group(2))

import llvmlite  # noqa: E401
import numba
import numpy as np
import pandas as pd
import pysam
import scipy
from scipy import stats

import dmrdetector


assert dmrdetector.__file__.startswith(
    "/opt/dmrdetector/src/"), dmrdetector.__file__
for mod in ("pip", "setuptools", "pytest"):
    assert importlib.util.find_spec(
        mod) is None, f"{mod} must not be in the image"

a = np.arange(10.0)
assert abs(stats.ttest_ind(a, a + 5).pvalue) < 0.05
assert pd.DataFrame({
    "g": list("aabb"),
    "x": [1, 2, 3, 4]
}).groupby("g").x.sum().tolist() == [3, 7]
assert np.linalg.inv(np.eye(3)).trace() == 3.0


# numba compiles and runs a kernel (cache off: the image is read-only)
@numba.njit(cache=False, nogil=True)
def _acc(x):
    s = 0.0
    for v in x:
        s += v * v
    return s


assert _acc(np.arange(5.0)) == 30.0

# kernels are compiled for a generic CPU, so the cache of the image is valid on
# any x86-64 host; the cache is compiled by env/warm_cache.sh, which leaves it
# out (with a WARNING in the build log) when a step fails, unless the build
# sets DMRDETECTOR_REQUIRE_CACHE=1
assert numba.config.CPU_NAME == "generic", numba.config.CPU_NAME

# glibc selects the FMA or the plain variant of its libm functions from the CPU
# (ifunc), and the two differ in the last bit of exp, sin and pow for 0.05% to
# 0.08% of the arguments. scipy.fft computes its twiddle factors and scipy.stats
# its distributions with them. The recipes set GLIBC_TUNABLES so that every host
# uses the plain variant. The expected values below are those of the plain
# variant of glibc 2.41 (the FMA variant gives 0x1.3ee6140bb23d3p+200 and
# 0x1.37a7f1c739a62p-1); they fail when the setting is lost or does not work.
if platform.machine() == "x86_64":
    assert os.environ.get("GLIBC_TUNABLES") == "glibc.cpu.hwcaps=-FMA,-FMA4", (
        os.environ.get("GLIBC_TUNABLES"))
    import math
    for fn, arg, want in ((math.exp, "0x1.15b2c1780da74p+7",
                           "0x1.3ee6140bb23d4p+200"),
                          (math.sin, "-0x1.073d468e2c0a9p+8",
                           "0x1.37a7f1c739a63p-1")):
        got = fn(float.fromhex(arg)).hex()
        assert got == float.fromhex(want).hex(), (
            f"libm {fn.__name__}({arg}) = {got}, expected the plain variant "
            f"{want}: GLIBC_TUNABLES has no effect")


def _write_samples(folder, n=6, sites=3000):
    """Writes n samples (BGZF 5-column BED and a one-sample VCF) on chr21."""
    import numpy as np
    import pysam
    rng = np.random.default_rng(1)
    pos = 5_010_000 + 2 * np.arange(sites)
    meth = np.where((pos // 20_000) % 2 == 0, 0.04, 0.85)  # blocks
    beds, vcfs = [], []
    for k in range(n):
        cov = rng.integers(10, 40, sites)
        nmod = rng.binomial(cov, meth)
        bed = os.path.join(folder, f"S{k}.bed")
        with open(bed, "w") as fh:
            for p, c, m in zip(pos, cov, nmod):
                fh.write(f"chr21\t{p}\t{p + 1}\t{c}\t{m}\n")
        pysam.tabix_compress(bed, bed + ".gz", force=True)
        pysam.tabix_index(bed + ".gz", preset="bed", force=True)
        vcf = os.path.join(folder, f"S{k}.vcf")
        with open(vcf, "w") as fh:
            fh.write(
                "##fileformat=VCFv4.2\n##contig=<ID=chr21>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS{k}\n"
            )
            for p in pos[100:300:10]:
                fh.write(f"chr21\t{p + 1}\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n")
        beds.append(bed + ".gz")
        vcfs.append(vcf)
    return beds, vcfs


def _run_kernels(folder):
    """Runs a small data set through the compiled steps of the library."""
    from dmrdetector.fileio import atlas_io
    from dmrdetector.fileio import bgzf
    from dmrdetector.steps import atlas
    from dmrdetector.steps import cpg_islands
    from dmrdetector.steps import preprocess

    beds, vcfs = _write_samples(folder)
    kept = []
    for k, (bed, vcf) in enumerate(zip(beds, vcfs)):
        out = os.path.join(folder, f"F{k}.bed.gz")
        res = preprocess.preprocess_sample(bed,
                                           vcf,
                                           5,
                                           out,
                                           log=lambda m: None,
                                           threads=1)
        assert res["rows_out"] > 0 and res["removed_by_het"] > 0, res
        kept.append(out)
    raw = os.path.join(folder, "raw_atlas.bed.gz")
    atlas.build_atlas(kept, raw, 5, 0.5, threads=2, log=lambda m: None)
    final = os.path.join(folder, "ref_atlas_150.bed.gz")
    atlas.smooth_atlas(raw, 150, final, False, 2, log=lambda m: None)
    assert bgzf.count_newlines(final, 1) > 1000
    islands = cpg_islands.detect(atlas_io.read_atlas(final),
                                 "autosomal",
                                 include_sex=False,
                                 log=lambda m: None)
    assert len(islands) > 0


def _cache_files(folder):
    """Returns {relative path: (size, mtime_ns)} of the cache files."""
    out = {}
    for top, _, names in os.walk(folder):
        for n in names:
            if n.endswith((".nbc", ".nbi")):
                st = os.stat(os.path.join(top, n))
                out[os.path.relpath(os.path.join(top, n),
                                    folder)] = (st.st_size, st.st_mtime_ns)
    return out


def _dispatchers():
    """Returns {id: (module name, function, hits, misses)} of the kernels."""
    from numba.core.dispatcher import Dispatcher
    found = {}
    for name, mod in list(sys.modules.items()):
        if not name.startswith("dmrdetector."):
            continue
        for obj in list(vars(mod).values()):
            if isinstance(obj, Dispatcher) and id(obj) not in found:
                st = obj.stats
                found[id(obj)] = (obj.py_func.__module__, obj.py_func.__name__,
                                  sum(st.cache_hits.values()),
                                  sum(st.cache_misses.values()))
    return found


def check_kernel_cache():
    """Asserts that the baked kernels load and that nothing is compiled."""
    seed = os.path.join(SEED_DIR, "SEED")
    if not os.path.isfile(seed):
        assert not REQUIRE_CACHE, (
            "DMRDETECTOR_REQUIRE_CACHE=1, but the image has no kernel cache "
            f"({seed} is missing): the build step env/warm_cache.sh failed")
        print("numba kernel cache: absent (not required; the first call "
              "compiles)")
        return
    modules = set()  # modules that have baked entries
    for top, _, names in os.walk(SEED_DIR):
        for n in names:
            if n.endswith(".nbi"):
                pkg = os.path.basename(top).rsplit("_", 1)[0]
                modules.add(f"dmrdetector.{pkg}.{n.split('.')[0]}")
    assert modules, "the kernel cache of the image has no index file"
    folder = os.environ["NUMBA_CACHE_DIR"]
    from dmrdetector.util import jit
    jit.set_default_cache_dir()  # copies the baked cache into `folder`
    seeded = _cache_files(folder)
    assert seeded, "the kernels were not copied from the cache of the image"
    data = os.path.join(WORK, "data")
    os.makedirs(data)
    _run_kernels(data)
    # a kernel compiled in a thread, a worker or a child would save a file
    assert _cache_files(folder) == seeded, "a kernel was compiled and saved"
    kernels = _dispatchers()
    misses = sorted(f"{m}.{f}" for m, f, _, miss in kernels.values() if miss)
    assert not misses, f"kernels compiled instead of loaded: {misses}"
    hits = {}
    for m, _, hit, _ in kernels.values():
        hits[m] = hits.get(m, 0) + hit
    assert sum(hits.values()) > 0, "no kernel was loaded from the cache"
    idle = sorted(m for m in modules if not hits.get(m))
    assert not idle, f"no kernel loaded from the cache in: {idle}"
    print("numba kernel cache:",
          open(seed).read().strip().replace("\n",
                                            "; "), "|", sum(hits.values()),
          "kernels loaded in", len(modules), "modules, 0 compiled")


check_kernel_cache()

# the BED files ship inside the package (src/dmrdetector/data, read with
# importlib.resources)
from dmrdetector import data  # noqa: E402


for f in (data.asset("par_hg38.bed"), data.asset("icr_gene_list.bed"),
          data.asset("target_region.hg38.bed")):
    assert os.path.isfile(f) and f.startswith("/opt/dmrdetector/src/"), f

# byte-compiled files ship in the image (PYTHONDONTWRITEBYTECODE is set at run
# time)
import importlib.util as iu  # noqa: E402


assert os.path.isfile(iu.cache_from_source(
    dmrdetector.__file__)), "dmrdetector is not byte-compiled"

# libdeflate: pinned Debian package, found where the code looks first
#   (sys.prefix/lib), round trip against zlib
so = os.path.join(sys.prefix, "lib", "libdeflate.so.0")
assert open(os.path.join(sys.prefix, "lib",
                         "libdeflate.version")).read().strip() == "1.23-2"
ld = ctypes.CDLL(so)
vp, sz = ctypes.c_void_p, ctypes.c_size_t
ld.libdeflate_alloc_compressor.restype = vp
ld.libdeflate_alloc_compressor.argtypes = [ctypes.c_int]
ld.libdeflate_alloc_decompressor.restype = vp
ld.libdeflate_deflate_compress.argtypes = [vp, vp, sz, vp, sz]
ld.libdeflate_deflate_compress.restype = sz
ld.libdeflate_deflate_decompress.argtypes = [vp, vp, sz, vp, sz, vp]
ld.libdeflate_deflate_decompress.restype = ctypes.c_int
raw = (b"chr1\t100\t101\t0.5\n" * 5000)
out = ctypes.create_string_buffer(len(raw))
c = ld.libdeflate_alloc_compressor(6)
n = ld.libdeflate_deflate_compress(c, raw, len(raw), out, len(raw))
assert 0 < n < len(raw)
assert zlib.decompress(out.raw[:n], -15) == raw  # zlib reads libdeflate output
back = ctypes.create_string_buffer(len(raw))
d = ld.libdeflate_alloc_decompressor()
assert ld.libdeflate_deflate_decompress(d, out, n, back, len(raw),
                                        None) == 0 and back.raw == raw
from dmrdetector.fileio import bgzf as _bgzf  # noqa: E402


assert _bgzf.libdeflate(), "dmrdetector does not find libdeflate"

print("python", platform.python_version(), "numpy", np.__version__, "scipy",
      scipy.__version__, "pandas", pd.__version__, "numba", numba.__version__,
      "llvmlite", llvmlite.__version__, "libdeflate 1.23-2, pysam",
      pysam.__version__, "dmrdetector", dmrdetector.__version__)
sys.exit(0)
