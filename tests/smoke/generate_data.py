#!/usr/bin/env python3
"""Writes the synthetic demo data set used by the README.

It reuses the data model of run_smoke.py: 8 controls (4 male, 4 female) and 2
cases (CASE01 female, CASE02 male) on chr22, chrX and chrY, with five CpG
islands and one ICR-like block on chr22 that are hypermethylated in the cases.
The data are tiny, so the folder also gets a settings.toml with relaxed
thresholds and a small synthetic Ensembl-style GTF.

Writes into OUTDIR (default ./demo):
  controls.csv            samplesheet of `build` (controls only)
  settings.toml           relaxed thresholds for the tiny data
  genes.gtf               synthetic genes at the islands (GRCh38 contig names)
  target_region.bed       one target region on the ICR-like block (the bundled
                          hg38 target regions lie outside the synthetic data)
  data/5col/              5-column BED (bgzip) and VCF of every sample
  data/modkit/            the same counts as modkit bedMethyl
  data/pbcpg/             the same counts as pb-CpG-tools

Needs numpy and pysam (both are in the image). Run it through the image:
  apptainer exec "$SIF" python tests/smoke/generate_data.py demo
"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_smoke  # noqa: E402  (needs the path above)

GENES = [  # (contig, start, end, symbol), 1-based, around the synthetic islands
    ("22", s - 500, s + run_smoke.ISLAND_LEN + 500, f"DEMO{i}")
    for i, s in enumerate(run_smoke.ISLANDS22, 1)
] + [("X", 29_999_000, 30_002_000, "DEMOX"),
     ("Y", 19_999_000, 20_002_000, "DEMOY")]

SETTINGS = """\
# Relaxed thresholds for the tiny demo data (8 controls, about 5,000 CpGs).
# The algorithm is unchanged. Real data use the defaults of
# `dmrdetector config init`.
[cpg_islands]
cpg_min_samples_per_site = 2
cpg_min_island_length = 100

[icr]
icr_min_sites = 5

[dmr]
dmr_cpg_min_sites = 10
dmr_cpg_diff_threshold = 0.3
dmr_fdr_threshold = 0.05
"""


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "demo").resolve()
    out.mkdir(parents=True, exist_ok=True)
    run_smoke.write_data(out)
    sheet = run_smoke.sheet(out / "controls.csv", out / "data")
    # paths relative to OUTDIR, so the commands run from that folder
    text = sheet.read_text().replace(f"{out}/", "")
    sheet.write_text(text)
    (out / "settings.toml").write_text(SETTINGS)
    with open(out / "genes.gtf", "w") as fh:
        for i, (c, s, e, sym) in enumerate(GENES, 1):
            fh.write(f"{c}\tdemo\tgene\t{s}\t{e}\t.\t+\t.\t"
                     f'gene_id "ENSGDEMO{i:05d}"; gene_name "{sym}"; '
                     'gene_biotype "protein_coding";\n')
    start = run_smoke.ICR_ISLAND
    end = start + run_smoke.ISLAND_LEN
    (out / "target_region.bed").write_text(
        "CHROM\tSTART\tEND\tREGION_NAME\tGENE_SYMBOL\tDISEASE_NAME\n"
        f"chr22\t{start}\t{end}\tDEMO:ICR\tDGCR6L\tDemo syndrome (hypo)\n")
    print(f"wrote the demo data set to {out}")


if __name__ == "__main__":
    main()
