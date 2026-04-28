# autosomal_sex_stratagy.md

# Chromosome Stratification Strategy in DMRdetector (Manuscript Version)

## Overview

The manuscript version of DMRdetector adopts a **ploidy-aware and X-inactivation-aware methylome framework**, rather than treating sex chromosomes as extensions of autosomal analysis. The genome is partitioned into biologically distinct compartments with compartment-specific reference atlases, CpG island definitions, and DMR detection thresholds.

This strategy recognizes that:

* Autosomes follow diploid methylation behavior.
* Male X chromosome is hemizygous.
* Female X chromosome is subject to X-inactivation and exhibits hemi-methylated regulatory structure.
* Sex chromosomes therefore require separate statistical treatment.

---

# 1. Genome Partitioning

The genome is partitioned into three major compartments:

## 1.1 Autosomal compartment

**Chromosomes:**

```text
chr1–chr22
```

### Reference atlas

A shared autosomal methylation reference atlas is constructed using all control samples.

### CpG island definition

CORE sites:

```math
REF_METH < 0.10
phi < 0.05
```

SHORE sites:

```math
REF_METH < 0.30
```

within 150 bp of a CORE site.

### Autosomal DMR thresholds

Regions retained as DMRs must satisfy:

* |Diff| > 40%
* > 30 CpG sites
* FDR < 0.01

This high-effect threshold targets diagnostic-scale methylation aberrations.

---

## 1.2 Sex chromosome compartment

Sex chromosomes are handled separately from autosomes.

Pseudoautosomal regions (PAR) are excluded from X-specific modeling.

```text
X analysis = non-PAR chrX only
```

---

# 2. Sex-Stratified X Chromosome Modeling

## 2.1 Male X chromosome (male_X)

### Region

```text
non-PAR chrX in males
```

### Reference atlas

A male-specific X chromosome reference atlas is constructed.

### CpG island model

Male X uses an autosomal-like hypomethylated island model:

```math
REF_METH < 0.10
```

### DMR thresholds

Same thresholds as autosomes:

* |Diff| > 40%
* > 30 CpGs
* FDR <0.01

### Biological rationale

Male X is hemizygous and does not exhibit X-inactivation mosaic methylation.

---

## 2.2 Female X chromosome (female_X)

### Region

```text
non-PAR chrX in females
```

### Reference atlas

A separate female-specific X chromosome atlas is constructed.

This is necessary because X-inactivation changes baseline methylation behavior.

---

## Female X CpG island definition

Unlike autosomes or male X:

```math
0.40 <= REF_METH <= 0.60
```

This hemi-methylated criterion models X-inactivation domains.

This replaces the canonical hypomethylated island definition.

---

## Female X DMR thresholds

Relaxed thresholds are used:

* |Diff| > 25%
* > 5 CpGs
* FDR <0.01

### Rationale

Because female X starts from a ~50% hemi-methylated baseline, biologically meaningful deviations are inherently smaller.

---

# 3. Y Chromosome

## Region

```text
chrY
```

Included only for male samples.

Treated as a male-specific compartment.

The manuscript places less emphasis on chrY-specific statistical modeling than on chrX, but conceptually it forms part of the sex-chromosome compartment.

---

# 4. Reference Atlases Produced

The manuscript framework implicitly generates three reference atlases:

```text
1. Autosomal atlas
2. Male chrX atlas
3. Female chrX atlas
```

(With chrY handled as a male-specific component.)

This is not a single pooled genome-wide reference.

---

# 5. Sample Routing Logic

## Male case samples

Run against:

```text
Autosomal atlas
Male X atlas
(optional chrY analysis)
```

---

## Female case samples

Run against:

```text
Autosomal atlas
Female X atlas
```

---

# 6. Summary of Statistical Regimes

| Compartment | Region       | Reference model | Island criterion       | DMR threshold |
| ----------- | ------------ | --------------- | ---------------------- | ------------- |
| Autosomal   | chr1-22      | Shared atlas    | REF_METH <10%          | 40%, >30 CpGs |
| Male X      | non-PAR chrX | Male-specific   | REF_METH <10%          | 40%, >30 CpGs |
| Female X    | non-PAR chrX | Female-specific | 40–60% hemi-methylated | 25%, >5 CpGs  |
| Y           | chrY         | Male-specific   | exploratory            | optional      |

---

# 7. Conceptual Framework

This strategy can be viewed as three statistical methylation regimes:

## Regime I

Autosomal diploid hypomethylated model

---

## Regime II

Male hemizygous X model

---

## Regime III

Female X-inactivation hemi-methylated model

---

Rather than treating sex chromosomes as special cases of autosomal analysis, the manuscript treats them as **distinct biological methylome domains requiring separate models.**

---

# 8. Reviewer-Style One-Sentence Summary

> The genome was partitioned into autosomal, male non-PAR X, female non-PAR X, and Y compartments, with female X modeled separately using hemi-methylated (40–60%) criteria to accommodate X-inactivation, and locus-specific DMR thresholds applied according to chromosomal biology.

---

# 9. Key Design Principle

This framework is best described as:

```text
Ploidy-aware + X-inactivation-aware hierarchical methylome modeling
```

rather than conventional sex-stratified analysis.

That distinction is central to the manuscript methodology.
