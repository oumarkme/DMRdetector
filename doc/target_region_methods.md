# Target regions for the DMRdetector v3 target-region report

`src/dmrdetector/data/assets/target_region.hg38.bed` lists 21 GRCh38 regions where a methylation change, or a lesion read out by methylation, is linked to a human disease. `src/dmrdetector/data/assets/target_region.sources.tsv` gives the coordinate basis, the direction source, the confidence, and notes for each row, plus five excluded candidates. Both files were built on 2026-10-06. Abbreviations in `sources.tsv`: AS, Angelman syndrome; BWS, Beckwith-Wiedemann syndrome; CA-DMR, clinically associated DMR; GOM and LOM, gain and loss of methylation; KOS, Kagami-Ogata syndrome; MBCS, Mulchandani-Bhoj-Conlin syndrome; PHP1B, pseudohypoparathyroidism type 1B; PWS, Prader-Willi syndrome; SRS, Silver-Russell syndrome; TNDM, transient neonatal diabetes mellitus; TS14, Temple syndrome; UPD, uniparental disomy.

## File format

The BED has a header line and six tab-separated columns: `CHROM`, `START` (0-based), `END`, `REGION_NAME`, `GENE_SYMBOL`, `DISEASE_NAME`. Genes or diseases that share a row are separated by `;`. Each disease carries one direction token in parentheses:

* `hyper` or `hypo`: the methylation change seen in the disease, relative to the normal hemi-methylated state.
* Any other text in parentheses marks a row without a methylation direction (AS-IC, UBE3A, IGF1R). Commas, never semicolons, are used inside the parentheses.

The direction is causal only for an epimutation. For uniparental disomy (UPD), deletion, or duplication, the same shift is the readout of the genetic lesion. The `NOTES` column says which applies. For DMPK the CTG expansion is the cause, and `hyper` describes the flanking methylation that comes with large expansions.

## Coordinates

We used three bases, in this order of preference.

1. Atlas-derived. In the manuscript atlas (69 controls, the archived `ref_atlas_150.bed.gz`), we took CpGs with N ≥ 5 samples and `REF_METH` 0.40-0.60. CpGs ≤ 150 bp apart were grouped, and clusters with ≥ 5 CpGs were kept. An internal 109-control atlas was scanned the same way as a check. Each cluster was compared with the published GRCh38 window of the same differentially methylated region (DMR) in Urakawa et al. (2025, Table S1). For a new row, `START`/`END` is the span of all manuscript-atlas clusters that overlap the published window.
2. Published hg38 window. Used for imprinting centers 1 and 2 at 11p15.5 (IC1 and IC2; Moriura et al., 2025), and for the AS imprinting center (AS-IC). For the AS-IC, the two primers that Beygo et al. (2020) used to amplify the whole AS-IC each match once in the GRCh38 primary chromosomes (exact match). The six AS-IC variants they list lie inside the resulting amplicon.
3. Analysis by the authors or annotation. The three DMPK regions come from an analysis script of the authors (positions read as 1-based inclusive; reading the positions as modkit 0-based starts selects the same CpGs). The UBE3A row is the Ensembl release 86 gene span.

The 12 rows carried over from an earlier internal version (DMRdetector2, not published) keep their coordinates. Its DMPK row was replaced. Two satellite clusters inside published windows were not merged into those rows: GRB10, chr7:50,781,399-50,781,644 (13 CpGs), and MEST, chr7:130,493,610-130,493,841 (10 CpGs).

## Regions

Coordinates are GRCh38, 1-based inclusive. CpGs is the number of manuscript-atlas sites in the region.

| Region | Coordinates | Disease (direction) | Basis | CpGs |
|---|---|---|---|---|
| PLAGL1:alt-TSS-DMR | chr6:144,007,046-144,008,785 | Transient neonatal diabetes mellitus (hypo) | atlas | 143 |
| GRB10:alt-TSS-DMR | chr7:50,781,905-50,783,354 | Silver-Russell syndrome (hyper) | atlas, v2 | 141 |
| PEG10:TSS-DMR | chr7:94,656,575-94,658,647 | Silver-Russell syndrome (hyper) | atlas, v2 | 92 |
| MEST:alt-TSS-DMR | chr7:130,490,281-130,493,449 | Silver-Russell syndrome (hyper) | atlas, v2 | 210 |
| H19/IGF2:IG-DMR | chr11:1,997,582-2,003,510 | Beckwith-Wiedemann (hyper); Silver-Russell (hypo) | published, v2 | 254 |
| IGF2:Ex9-DMR | chr11:2,132,761-2,133,881 | Beckwith-Wiedemann (hyper); Silver-Russell (hypo) | atlas, v2 | 63 |
| IGF2:alt-TSS-DMR | chr11:2,146,415-2,148,809 | Beckwith-Wiedemann (hyper); Silver-Russell (hypo) | atlas, v2 | 59 |
| KCNQ1OT1:TSS-DMR | chr11:2,698,718-2,701,029 | Beckwith-Wiedemann (hypo) | published, v2 | 193 |
| MEG3/DLK1:IG-DMR | chr14:100,809,368-100,811,319 | Temple (hypo); Kagami-Ogata (hyper) | atlas, v2 | 53 |
| MEG3:TSS-DMR | chr14:100,824,147-100,827,640 | Temple (hypo); Kagami-Ogata (hyper) | atlas, v2 | 191 |
| SNRPN:AS-IC | chr15:24,920,260-24,921,630 | Angelman (maternal AS-IC deletion, not a DMR) | published primers | 19 |
| SNURF:TSS-DMR | chr15:24,954,857-24,956,828 | Prader-Willi (hyper); Angelman (hypo) | atlas, v2 | 113 |
| UBE3A:gene | chr15:25,333,728-25,439,056 | Angelman (locus, methylation not diagnostic) | Ensembl 86 | 844 |
| IGF1R:Int2-DMR | chr15:98,865,177-98,866,881 | IGF1R-related growth restriction (direction not established) | atlas, v2 | 67 |
| DMPK:Region-C | chr19:45,769,405-45,770,205 | Myotonic dystrophy type 1 (hyper) | authors | 66 |
| DMPK:Region-B | chr19:45,770,262-45,770,662 | Myotonic dystrophy type 1 (hyper) | authors | 41 |
| DMPK:Region-A | chr19:45,770,752-45,771,042 | Myotonic dystrophy type 1 (hyper) | authors | 28 |
| GNAS-NESP:TSS-DMR | chr20:58,839,107-58,843,373 | Pseudohypoparathyroidism type 1B (hyper); Mulchandani-Bhoj-Conlin (hypo) | atlas | 255 |
| GNAS-AS1:TSS-DMR | chr20:58,850,397-58,852,286 | Pseudohypoparathyroidism type 1B (hypo); Mulchandani-Bhoj-Conlin (hyper) | atlas | 93 |
| GNAS-XL:Ex1-DMR | chr20:58,853,850-58,856,416 | Pseudohypoparathyroidism type 1B (hypo); Mulchandani-Bhoj-Conlin (hyper) | atlas | 201 |
| GNAS-A/B:TSS-DMR | chr20:58,888,210-58,890,027 | Pseudohypoparathyroidism type 1B (hypo); Mulchandani-Bhoj-Conlin (hyper) | atlas | 187 |

Five candidates were excluded. The DMPK lowly methylated region of the earlier version was replaced. The NDN, MAGEL2 and MKRN3 TSS-DMRs are non-clinical DMRs (Urakawa et al., 2025), and no disease caused by methylation change there was found. The 15q11-q13 deletion breakpoints BP1-BP3 have no hg38 coordinates in the sources reviewed and are copy-number targets.

## Changes from the earlier internal version (DMRdetector2)

* MEG3/DLK1:IG-DMR and MEG3:TSS-DMR were labeled Silver-Russell syndrome. Both rows now carry Temple syndrome (hypo) and Kagami-Ogata syndrome (hyper) (Ogata & Kagami, 2016; Urakawa et al., 2025).
* SNURF:TSS-DMR gained Angelman syndrome (hypo).
* IGF1R:Int2-DMR is no longer labeled Silver-Russell syndrome. Its README called it a clinically associated DMR. Urakawa et al. (2025, Table 2) list it, and PEG10:TSS-DMR, as non-clinical DMRs.
* New rows: PLAGL1, the AS-IC, UBE3A, three DMPK regions, and four GNAS DMRs.

## Limitations

* Atlas-derived boundaries come from blood of 69 controls. These boundaries were not copied from a paper, and they differ from the published windows by up to about 700 bp (IGF2:alt-TSS-DMR).
* The AS-IC and UBE3A rows cannot give a methylation diagnosis. The AS-IC is about 79% methylated in the atlas and is not a DMR (Matsubara et al., 2019); a maternal AS-IC deletion is detected by copy number and by hypomethylation at SNURF:TSS-DMR. UBE3A sequence variants (5-10% of Angelman syndrome) need sequencing (Beygo et al., 2019).
* IGF1R:Int2-DMR is kept by decision of the authors; no methylation-defined disease is established there.
* The target report still assigns DIRECTION (hyper, hypo, normal) to the three rows without a disease direction. UBE3A:gene (844 CpGs) and SNRPN:AS-IC (19 CpGs) overlap genes in `src/dmrdetector/data/assets/icr_gene_list.bed` and get the ICR cut-off (difference 0.25); IGF1R:Int2-DMR gets the non-ICR cut-off. Ignore DIRECTION on these rows.
* At DMPK the merged signal of a methylated expanded allele is at most about 0.5, and reads through a long expansion may not align. A DIRECTION call at the 0.40 cut-off will therefore be rare, even when the flanks are methylated.
* DMPK and IGF1R do not overlap the ICR panel, so the target report applies the non-ICR cut-off (difference 0.40). The report has no minimum CpG count: a region needs 2 CpGs for the t-test and otherwise has `NA` in `P_FDR` and `DIRECTION`. DMPK:Region-A (28 atlas CpGs) gets a DIRECTION value like the other regions.
* The evidence for IGF2 DMR2 in Silver-Russell syndrome rests on a small series (Higashimoto et al., 2021). IGF2 DMR0 loses methylation with age in healthy people (Sun et al., 2019).
* Which DMPK region holds the CTCF1 or CTCF2 site of Barbé et al. (2017) was not verified.

## References

Barbé, L., Lanni, S., López-Castel, A., Franck, S., Spits, C., Keymolen, K., Seneca, S., Tomé, S., Miron, I., Letourneau, J., Liang, M., Choufani, S., Weksberg, R., Wilson, M. D., Sedlacek, Z., Gagnon, C., Musova, Z., Chitayat, D., Shannon, P., ... Pearson, C. E. (2017). CpG methylation, a parent-of-origin effect for maternal-biased transmission of congenital myotonic dystrophy. *American Journal of Human Genetics, 100*(3), 488-505. https://doi.org/10.1016/j.ajhg.2017.01.033

Beygo, J., Buiting, K., Ramsden, S. C., Ellis, R., Clayton-Smith, J., & Kanber, D. (2019). Update of the EMQN/ACGS best practice guidelines for molecular analysis of Prader-Willi and Angelman syndromes. *European Journal of Human Genetics, 27*(9), 1326-1340. https://doi.org/10.1038/s41431-019-0435-0

Beygo, J., Grosser, C., Kaya, S., Mertel, C., Buiting, K., & Horsthemke, B. (2020). Common genetic variation in the Angelman syndrome imprinting centre affects the imprinting of chromosome 15. *European Journal of Human Genetics, 28*(6), 835-839. https://doi.org/10.1038/s41431-020-0595-y

Eggermann, K., Bliek, J., Brioude, F., Algar, E., Buiting, K., Russo, S., Tümer, Z., Monk, D., Moore, G., Antoniadi, T., Macdonald, F., Netchine, I., Lombardi, P., Soellner, L., Begemann, M., Prawitt, D., Maher, E. R., Mannens, M., Riccio, A., ... Eggermann, T. (2016). EMQN best practice guidelines for the molecular genetic testing and reporting of chromosome 11p15 imprinting disorders: Silver-Russell and Beckwith-Wiedemann syndrome. *European Journal of Human Genetics, 24*(10), 1377-1387. https://doi.org/10.1038/ejhg.2016.45

Han, Y., Jang, J. H., & Chang, H. (2026). Targeted long-read sequencing for high-resolution repeat profiling in myotonic dystrophy type 1. *Experimental & Molecular Medicine, 58*(4), 1203-1215. https://doi.org/10.1038/s12276-026-01683-6

Hannula-Jouppi, K., Muurinen, M., Lipsanen-Nyman, M., Reinius, L. E., Ezer, S., Greco, D., & Kere, J. (2014). Differentially methylated regions in maternal and paternal uniparental disomy for chromosome 7. *Epigenetics, 9*(3), 351-365. https://doi.org/10.4161/epi.27160

Higashimoto, K., Watanabe, H., Tanoue, Y., Tonoki, H., Tokutomi, T., Hara, S., Yatsuki, H., & Soejima, H. (2021). Hypomethylation of a centromeric block of ICR1 is sufficient to cause Silver-Russell syndrome. *Journal of Medical Genetics, 58*(6), 422-425. https://doi.org/10.1136/jmedgenet-2020-106907

Maeda, T., Higashimoto, K., Jozaki, K., Yatsuki, H., Nakabayashi, K., Makita, Y., Tonoki, H., Okamoto, N., Takada, F., Ohashi, H., Migita, M., Kosaki, R., Matsubara, K., Ogata, T., Matsuo, M., Hamasaki, Y., Ohtsuka, Y., Nishioka, K., Joh, K., ... Soejima, H. (2014). Comprehensive and quantitative multilocus methylation analysis reveals the susceptibility of specific imprinted differentially methylated regions to aberrant methylation in Beckwith-Wiedemann syndrome with epimutations. *Genetics in Medicine, 16*(12), 903-912. https://doi.org/10.1038/gim.2014.46

Mantovani, G., Bastepe, M., Monk, D., de Sanctis, L., Thiele, S., Usardi, A., Ahmed, S. F., Bufo, R., Choplin, T., De Filippo, G., Devernois, G., Eggermann, T., Elli, F. M., Freson, K., García Ramirez, A., Germain-Lee, E. L., Groussin, L., Hamdy, N., Hanna, P., ... Linglart, A. (2018). Diagnosis and management of pseudohypoparathyroidism and related disorders: First international consensus statement. *Nature Reviews Endocrinology, 14*(8), 476-500. https://doi.org/10.1038/s41574-018-0042-0

Matsubara, K., Itoh, M., Shimizu, K., Saito, S., Enomoto, K., Nakabayashi, K., Hata, K., Kurosawa, K., Ogata, T., Fukami, M., & Kagami, M. (2019). Exploring the unique function of imprinting control centers in the PWS/AS-responsible region: Finding from array-based methylation analysis in cases with variously sized microdeletions. *Clinical Epigenetics, 11*, 36. https://doi.org/10.1186/s13148-019-0633-1

McCullough, M. E., Letourneau-Freiberg, L. R., Bowden, T. L., Kandasamy, B., Ray, A., Wroblewski, K., Del Gaudio, D., Mackay, D. J. G., Philipson, L. H., & Greeley, S. A. W. (2024). Clinical characteristics and remission monitoring of 6q24-related transient neonatal diabetes. *Pediatric Diabetes, 2024*, 3624339. https://doi.org/10.1155/pedi/3624339

Moriura, Y., Nishio, Y., Ichimura, S., Noda, H., Tanahashi, Y., Yamamoto, H., Nakazawa, Y., Oso, T., Sato, Y., Takenouchi, T., Saitoh, S., Muramatsu, Y., & Ogi, T. (2025). De novo *CDKN1C* variant in Beckwith-Wiedermann spectrum with atypical complications. *Human Genome Variation, 12*(1), 9. https://doi.org/10.1038/s41439-025-00316-0

Mulchandani, S., Bhoj, E. J., Luo, M., Powell-Hamilton, N., Jenny, K., Gripp, K. W., Elbracht, M., Eggermann, T., Turner, C. L., Temple, I. K., Mackay, D. J., Dubbs, H., Stevenson, D. A., Slattery, L., Zackai, E. H., Spinner, N. B., Krantz, I. D., & Conlin, L. K. (2016). Maternal uniparental disomy of chromosome 20: A novel imprinting disorder of growth failure. *Genetics in Medicine, 18*(4), 309-315. https://doi.org/10.1038/gim.2015.103

Murrell, A., Ito, Y., Verde, G., Huddleston, J., Woodfine, K., Silengo, M. C., Spreafico, F., Perotti, D., De Crescenzo, A., Sparago, A., Cerrato, F., & Riccio, A. (2008). Distinct methylation changes at the IGF2-H19 locus in congenital growth disorders and cancer. *PLoS ONE, 3*(3), e1849. https://doi.org/10.1371/journal.pone.0001849

Ogata, T., & Kagami, M. (2016). Kagami-Ogata syndrome: A clinically recognizable upd(14)pat and related disorder affecting the chromosome 14q32.2 imprinted region. *Journal of Human Genetics, 61*(2), 87-94. https://doi.org/10.1038/jhg.2015.113

Prasasya, R., Grotheer, K. V., Siracusa, L. D., & Bartolomei, M. S. (2020). Temple syndrome and Kagami-Ogata syndrome: Clinical presentations, genotypes, models and mechanisms. *Human Molecular Genetics, 29*(R1), R107-R116. https://doi.org/10.1093/hmg/ddaa133

Sun, F., Higashimoto, K., Awaji, A., Ohishi, K., Nishizaki, N., Tanoue, Y., Aoki, S., Watanabe, H., Yatsuki, H., & Soejima, H. (2019). The extent of DNA methylation anticipation due to a genetic defect in ICR1 in Beckwith-Wiedemann syndrome. *Journal of Human Genetics, 64*(9), 937-943. https://doi.org/10.1038/s10038-019-0634-0

Tannorella, P., Minervino, D., Guzzetti, S., Vimercati, A., Calzari, L., Patti, G., Maghnie, M., Allegri, A. E. M., Milani, D., Scuvera, G., Mariani, M., Modena, P., Selicorni, A., Larizza, L., & Russo, S. (2021). Maternal uniparental disomy of chromosome 20 (UPD(20)mat) as differential diagnosis of Silver Russell syndrome: Identification of three new cases. *Genes, 12*(4), 588. https://doi.org/10.3390/genes12040588

Urakawa, T., Hattori, A., Ogiwara, Y., Masubuchi, H., Igarashi, M., Nakamura, S., Hara-Isono, K., Ishiwata, K., Ogata-Kawata, H., Kamura, H., Kuroki, Y., Nakabayashi, K., Fukami, M., & Kagami, M. (2025). A comprehensive long-read sequencing system to assess DNA methylation at differentially methylated regions and imprinting-disorder-related genes. *Genome Medicine, 17*, 144. https://doi.org/10.1186/s13073-025-01559-w
