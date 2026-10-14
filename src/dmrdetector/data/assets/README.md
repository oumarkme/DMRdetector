# Bundled assets

Files that `build`, `run` and the stage commands read when no option names another file. They ship inside the Python package (`src/dmrdetector/data/assets/`) and are read through `importlib.resources`. Provenance is what the repository records (commit dates are from `git log`); anything else is marked "not recorded".

| File | Content | Read by | Provenance |
|---|---|---|---|
| `icr_gene_list.bed` | Imprinting gene candidate intervals, header `CHROM START END GENE`; 124 rows, 122 gene symbols (`SNORD116` three times) | `build` (ICR panel), `run` (ICR gene override in the annotation, target report), `icr`, `annotate` | List of imprinted genes from geneimprint.com, each gene region extended by 5 kb (`README.md`, ICR section). It is the 124-row list of commit 2bc9d61 (v2.2, 2026-04-29), with md5 1b9166e53b15751eeee7efa2c63b2dcb. Source annotation release: not recorded. |
| `par_hg38.bed` | Pseudoautosomal regions of GRCh38, PAR1 and PAR2 on X and Y, each with and without the `chr` prefix, no header | Preprocessing (sex split) | Added in commit a00657d (2026-04-27). Source of the coordinates: not recorded. |
| `target_region.hg38.bed` | 21 disease-linked methylation regions for `run --target-region`, with a direction token in `DISEASE_NAME` | `run --target-region`, `target` | Built 2026-10-06 (commit 9254a34); coordinates are atlas-derived or taken from the published windows named in `target_region.sources.tsv`. |
| `target_region.sources.tsv` | Coordinate basis, direction source, confidence and notes per region, and five excluded candidates | Documentation only | Same date; references in `doc/target_region_methods.md`. |

## ICR gene list

`icr_gene_list.bed` reproduces the manuscript panel of imprinting control regions (ICRs). With it, v2.2 (git tag `v2.2-manuscript`, unmodified) gives the 446-line `icr_region.bed` of the manuscript run (445 ICRs) and the six annotated tables of differentially methylated regions (DMRs) of the cases line by line. v3 gives the same regions, genes and directions; `P_VALUE` and `P_FDR` agree within 6.7e-13 relative (`doc/validation.md`).

Two properties of the list matter for anyone who edits it (see `doc/ALGORITHM.md`):

* The extension of 5 kb on both sides is part of the file. The code uses START and END as given and adds nothing.
* Of the 124 rows, 40 equal the Ensembl release 86 gene coordinates plus 5 kb each side. The other 84 do not (for example *TP73* differs by 4 bp at the start and *ANO1* by about 97 kb), and the annotation they come from is not recorded.

Until 2026-10-06 the file held 121 rows (commit 7e83486 had removed *MESTIT1*, *PARD6G-AS1* and *GNAS-AS1*); that list gives 425 ICRs on the manuscript atlas. The v2.2 Nextflow workflow (git tag `v2.2-manuscript`, removed from this tree) read `assets/icr_gene_list.bed` of its own tree by default.
