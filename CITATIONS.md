# Citations

## DMRdetector

> Ou J-H, Tsai M-JM, Hsu R-H, Chen H-A, Wu C-S, Chao M-Y, Hsu C, Chien Y-H, Lee N-C, Hwu W-L. Detecting differentially methylated regions at CpG islands and imprinted loci with an ONT methylation atlas. (Manuscript in preparation, 2026.) Code: https://github.com/oumarkme/DMRdetector

## Reference atlas

The manuscript atlas files (a Zenodo dataset, not an archive of the code): https://doi.org/10.5281/zenodo.18978867

## Software

- **NumPy**
  > Harris CR, et al. Array programming with NumPy. *Nature* 585, 357-362 (2020). https://doi.org/10.1038/s41586-020-2649-2

- **SciPy** (FFT, t distribution)
  > Virtanen P, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17, 261-272 (2020). https://doi.org/10.1038/s41592-019-0686-2

- **pandas**
  > McKinney W. Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56-61 (2010). https://doi.org/10.25080/Majora-92bf1922-00a

- **Numba**
  > Lam SK, Pitrou A, Seibert S. Numba: a LLVM-based Python JIT compiler. *Proceedings of the Second Workshop on the LLVM Compiler Infrastructure in HPC*, 1-6 (2015). https://doi.org/10.1145/2833157.2833162

- **pysam, htslib and bcftools** (VCF and tabix access; `pysam` wraps htslib)
  > Bonfield JK, et al. HTSlib: C library for reading/writing high-throughput sequencing data. *GigaScience* 10(2), giab007 (2021). https://doi.org/10.1093/gigascience/giab007
  > Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience* 10(2), giab008 (2021). https://doi.org/10.1093/gigascience/giab008

- **tabix**
  > Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. *Bioinformatics* 27(5), 718-719 (2011). https://doi.org/10.1093/bioinformatics/btq671

- **libdeflate** (block decompression and compression)
  > Biggers E. libdeflate. https://github.com/ebiggers/libdeflate

- **Apptainer / Singularity** (container runtime)
  > Kurtzer GM, Sochat V, Bauer MW. Singularity: scientific containers for mobility of compute. *PLoS ONE* 12(5), e0177459 (2017). https://doi.org/10.1371/journal.pone.0177459

- **Ensembl** (GRCh38 gene models; the manuscript used release 86)
  > Martin FJ, et al. Ensembl 2023. *Nucleic Acids Research* 51(D1), D933-D941 (2023). https://doi.org/10.1093/nar/gkac958

## Data sources

- **Imprinted gene list.** geneimprint.com (https://www.geneimprint.com), as used in the manuscript; see `src/dmrdetector/data/assets/README.md`.
- **Target regions.** The references for each region and direction are in `doc/target_region_methods.md` and `src/dmrdetector/data/assets/target_region.sources.tsv`.

## Methods

The statistical methods and their sources are listed in the manuscript and summarized in `doc/ALGORITHM.md`: Beta-Binomial dispersion (Feng et al., 2014), dispersion shrinkage (Love et al., 2014), single-case methylation testing (Rezwan et al., 2015; Grolaux et al., 2022), and the CORE/SHORE classification (Burger et al., 2013).
