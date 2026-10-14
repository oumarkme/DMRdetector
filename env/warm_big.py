"""Writes one large methylation BED per format for env/warm_cache.sh.

Real inputs are many megabytes, and the reader passes pieces of 4 MiB to the
numba kernels. A line that straddles two pieces is joined as `bytes`, so the
kernel then gets a read-only array, which is another compiled signature than
the writable array of the small demo files. The cache of the image must hold
both, so the demo needs files of more than one piece.

usage: warm_big.py SMOKE_DIR OUT_DIR
  SMOKE_DIR  folder with run_smoke.py (tests/smoke), for the line formats
  OUT_DIR    receives BIG.5col.bed.gz, BIG.modkit.bed.gz and BIG.pbcpg.bed.gz
             (bgzip), 200,000 CpGs each on chr22
"""

from pathlib import Path
import sys

import pysam

sys.path.insert(0, sys.argv[1])
import run_smoke  # noqa: E402  (needs the path above)

N_SITES = 200_000
FIRST = 20_000_000  # chr22, inside the contig; every second base is a CpG


def main() -> None:
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    for fmt, (line, header) in run_smoke.FORMATS.items():
        plain = out / f"BIG.{fmt}.bed"
        with open(plain, "w") as fh:
            fh.write(header)
            for i in range(N_SITES):
                pos = FIRST + 2 * i
                fh.write(line("chr22", pos, pos + 1, 25, 20))
        pysam.tabix_compress(str(plain), f"{plain}.gz", force=True)
        plain.unlink()


if __name__ == "__main__":
    main()
