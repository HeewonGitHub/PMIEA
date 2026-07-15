# PMIEA

**PMIEA: Phenotype-Predictive Molecular Interplays Enrichment Analysis using Graph-Guided Random Forest**

---

## Overview

This repository provides the R implementation and toy dataset for **Phenotype-Predictive Molecular Interplays Enrichment Analysis (PMIEA)**.

PMIEA is a gene network enrichment framework designed to identify biological pathways enriched with molecular interplays that contribute to phenotype discrimination. The method combines:

1. Graph-guided random forest construction
2. Gini impurity-based molecular interplay importance
3. Differential network topology
4. Edge-based pathway enrichment analysis
5. Gene-level permutation testing

The repository contains an example analysis comparing tumor and normal gene expression profiles.

---

## Repository Structure

```text
PMIEA/
├── README.md
├── sampleR_PMIEA.r
└── ToyDATA_PMIEA/
    ├── EXP_Tumor.csv
    ├── EXP_Normal.csv
    ├── NP_Tumor.csv
    ├── NP_Normal.csv
    └── Pathway_Genes.csv
```

---

## Software Requirements

The implementation is written in **R**.

Required packages:

```r
install.packages("HDeconometrics")
install.packages("Hmisc")
install.packages("igraph")
install.packages("bayesmeta")
```

Load the required packages:

```r
library(HDeconometrics)
library(Hmisc)
library(igraph)
library(bayesmeta)
```

---

## Input Data

The toy dataset is stored in the `ToyDATA_PMIEA` directory.

### Gene Expression Data

```text
ToyDATA_PMIEA/EXP_Tumor.csv
ToyDATA_PMIEA/EXP_Normal.csv
```

Rows represent samples and columns represent genes.

- `EXP_Tumor.csv`: gene expression matrix for the tumor group
- `EXP_Normal.csv`: gene expression matrix for the normal group

The two expression matrices must contain the same gene columns in the same order.

### Gene Network Data

```text
ToyDATA_PMIEA/NP_Tumor.csv
ToyDATA_PMIEA/NP_Normal.csv
```

The first two columns contain the two endpoint genes of each molecular interplay.

- `NP_Tumor.csv`: gene network inferred from tumor samples
- `NP_Normal.csv`: gene network inferred from normal samples

### Pathway Gene Set

```text
ToyDATA_PMIEA/Pathway_Genes.csv
```

This file contains the genes belonging to the pathway evaluated by PMIEA.

---

## PMIEA Workflow

```text
Tumor and normal gene expression data
                ↓
Tumor and normal gene network construction
                ↓
Graph-guided random forest
                ↓
Gini impurity-based molecular interplay importance
                ↓
Differential edge betweenness and hubness
                ↓
Topology-aware molecular interplay predictive score
                ↓
Edge-based pathway enrichment analysis
                ↓
Gene-level permutation test
                ↓
Permutation p-value
```

---

## Step 1: Prepare Phenotype Data

The tumor and normal expression matrices are combined, and phenotype labels are assigned to the samples.

```r
EXP_PMIEA <- rbind(EXP_TM, EXP_NM)

Y_PMIEA <- factor(
  c(
    rep("Tumor", nrow(EXP_TM)),
    rep("Normal", nrow(EXP_NM))
  )
)
```

---

## Step 2: Construct the Prior Gene Network

The tumor network is used as the prior graph for graph-guided random forest construction.

```r
g_prior <- graph_from_data_frame(
  EDGE_PRIOR,
  directed = FALSE,
  vertices = common_genes
)

g_prior <- simplify(g_prior)
```

---

## Step 3: Graph-Guided Random Forest

For each graph-guided decision tree:

1. A seed gene is randomly selected.
2. A local subnetwork is sampled using breadth-first search.
3. The root split is selected from genes in the sampled subnetwork.
4. Subsequent split genes are restricted to genes connected to the parent split gene.
5. Gini impurity reductions are accumulated for genes and molecular interplays.

The default parameters in the example script are:

```r
Rtree <- 100
sub_size <- 50
max_depth <- 5
min_node_size <- 10
```

---

## Step 4: Molecular Interplay Importance

For each molecular interplay, PMIEA aggregates the Gini impurity reductions of its connected genes across graph-guided decision trees.

The resulting molecular interplay importance values are min-max normalized to the interval from 0 to 1.

---

## Step 5: Topology-Aware Predictive Score

PMIEA incorporates differential network topology using:

- Differential edge betweenness
- Differential edge hubness
- Phenotype-predictive molecular interplay importance

The topology-aware molecular interplay predictive score is computed using a logistic transformation:

```r
IuvSCORE$W_uv <- 1 / (
  1 + exp(
    -(
      IuvSCORE$dB_uv_tilde +
      IuvSCORE$dH_uv_tilde +
      IuvSCORE$Iuv_norm
    )
  )
)
```

---

## Step 6: Pathway Enrichment Analysis

Molecular interplays are ranked according to their topology-aware predictive scores.

A pathway edge is defined as an edge for which both endpoint genes belong to the pathway gene set. The enrichment score measures whether pathway edges are preferentially concentrated near the top of the ranked molecular interplay list.

---

## Step 7: Statistical Significance

Statistical significance is evaluated using gene-level permutation.

The example analysis uses 200 permutations:

```r
NOpm <- 201
```

The first iteration corresponds to the observed pathway, and the remaining 200 iterations correspond to randomly sampled gene sets of the same size.

The empirical permutation p-value is printed at the end of the analysis:

```r
print(Pvalue)
```

---

## Running the Example

Place the following files in the same repository structure:

```text
sampleR_PMIEA.r
ToyDATA_PMIEA/
```

Run the analysis in R:

```r
source("sampleR_PMIEA.r")
```

Alternatively, run it from a terminal:

```bash
Rscript sampleR_PMIEA.r
```

The file paths in the example script use Windows-style separators:

```r
"ToyDATA_PMIEA\\NP_Tumor.csv"
```

For macOS or Linux, replace them with forward slashes:

```r
"ToyDATA_PMIEA/NP_Tumor.csv"
```

---

## Output

The main objects generated by the script are:

| Object | Description |
|---|---|
| `gene_importance_result` | Phenotype-predictive importance of individual genes |
| `IuvSCORE` | Molecular interplay importance and topology-aware predictive scores |
| `SC_T` | Observed and permutation-based enrichment scores |
| `nmSC_T` | Normalized enrichment scores |
| `Pvalue` | Empirical permutation p-value for pathway enrichment |

The final p-value is displayed in the R console.

---

## Reproducibility

The repository provides the R script and toy dataset required to reproduce the example PMIEA analysis.

Because graph-guided tree construction, bootstrap sampling, and pathway permutation involve random sampling, users may add the following command before running the analysis:

```r
set.seed(1234)
```

Researchers can apply PMIEA to their own datasets by replacing the toy expression matrices, gene networks, and pathway gene set while preserving the required input formats.

---

## Citation

If you use this repository, please cite:

> Park H, Imoto S.  
> **PMIEA: Phenotype-Predictive Molecular Interplays Enrichment Analysis using Graph-Guided Random Forest.**  
> Manuscript in preparation.

The citation information should be updated after publication.

---

## Contact

For questions regarding the PMIEA method or implementation, please contact:

**Heewon Park**  
Sungshin Women's University  
Email: heewonn.park@gmail.com
