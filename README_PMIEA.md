# PMIEA

**PMIEA: Phenotype-Predictive Molecular Interplays Enrichment Analysis using Graph-Guided Random Forest**

---

## Overview

This repository provides the R implementation and toy dataset for **Phenotype-Predictive Molecular Interplays Enrichment Analysis (PMIEA)**.

PMIEA identifies biological pathways enriched with molecular interplays that contribute to phenotype discrimination by integrating graph-guided random forests, Gini impurity-based molecular interplay importance, differential network topology, and edge-based enrichment analysis.

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

## Software Requirements

The analysis is implemented in **R**.

```r
install.packages("HDeconometrics")
install.packages("Hmisc")
install.packages("igraph")
install.packages("bayesmeta")
```

---

## Running the Example

Run the analysis from R:

```r
source("sampleR_PMIEA.r")
```

or from a terminal:

```bash
Rscript sampleR_PMIEA.r
```

The current script uses Windows-style paths. For macOS or Linux, replace `\\` with `/`.

---


# Complete PMIEA Analysis Code


## Required R Packages

```r
library(HDeconometrics)
library(Hmisc)
library(igraph)
library(bayesmeta)
```

---


# Step 0: Load Input Data

```r
# Network of Tumor group
NP_TM<-read.table("ToyDATA_PMIEA\\NP_Tumor.csv",sep=",",header=TRUE) 
# Network of Normal group
NP_NM<-read.table("ToyDATA_PMIEA\\NP_Normal.csv",sep=",",header=TRUE)

# Expression levels of Tumor group
EXP_TM<-read.table("ToyDATA_PMIEA\\EXP_Tumor.csv",sep=",")
# Expression levels of Normal group
EXP_NM<-read.table("ToyDATA_PMIEA\\EXP_Normal.csv",sep=",")

# Pathway Genes
Observed_Pathway_Genes <- read.table("ToyDATA_PMIEA\\Pathway_Genes.csv",sep=",")[,1]
```

---


# Step 1: Prepare Phenotype Data

```r
EXP_PMIEA <- rbind(EXP_TM, EXP_NM)

Y_PMIEA <- factor(
  c(rep("Tumor", nrow(EXP_TM)),
    rep("Normal", nrow(EXP_NM)))
)
common_genes <- colnames(EXP_PMIEA)
```

---


# Step 2: Construct the Prior Gene Network

```r
EDGE_PRIOR <- NP_TM[, 1:2, drop = FALSE]
EDGE_PRIOR <- as.data.frame(EDGE_PRIOR)
colnames(EDGE_PRIOR) <- c("u", "v")

EDGE_PRIOR$u <- as.character(EDGE_PRIOR$u)
EDGE_PRIOR$v <- as.character(EDGE_PRIOR$v)

EDGE_PRIOR <- EDGE_PRIOR[
  EDGE_PRIOR$u %in% common_genes &
    EDGE_PRIOR$v %in% common_genes,
]

g_prior <- graph_from_data_frame(
  EDGE_PRIOR,
  directed = FALSE,
  vertices = common_genes
)

g_prior <- simplify(g_prior)
```

---


# Step 3: Define Utility Functions

```r
gini_impurity <- function(y) {
  tab <- table(y)
  p <- tab / sum(tab)
  return(1 - sum(p^2))
}

best_split <- function(X, y, candidate_genes, min_node_size = 5) {
  
  parent_gini <- gini_impurity(y)
  n <- length(y)
  
  best_gene <- NULL
  best_cut <- NULL
  best_delta <- 0
  
  for (gene in candidate_genes) {
    
    x <- X[, gene]
    
    if (length(unique(x)) <= 2) next
    
    cuts <- unique(as.numeric(quantile(
      x,
      probs = seq(0.1, 0.9, by = 0.1),
      na.rm = TRUE
    )))
    
    for (cut in cuts) {
      
      left_id <- which(x <= cut)
      right_id <- which(x > cut)
      
      if (length(left_id) < min_node_size ||
          length(right_id) < min_node_size) next
      
      g_left <- gini_impurity(y[left_id])
      g_right <- gini_impurity(y[right_id])
      
      child_gini <-
        length(left_id) / n * g_left +
        length(right_id) / n * g_right
      
      delta <- parent_gini - child_gini
      
      if (delta > best_delta) {
        best_delta <- delta
        best_gene <- gene
        best_cut <- cut
      }
    }
  }
  
  return(list(
    gene = best_gene,
    cut = best_cut,
    delta = best_delta
  ))
}

get_bfs_subnetwork <- function(g, seed_gene, m = 50) {
  
  bfs_res <- bfs(
    g,
    root = seed_gene,
    order = TRUE,
    unreachable = FALSE
  )
  
  genes <- V(g)$name[bfs_res$order]
  genes <- genes[!is.na(genes)]
  genes <- genes[genes %in% V(g)$name]
  
  if (length(genes) > m) {
    genes <- genes[1:m]
  }
  
  return(genes)
}

build_graph_guided_tree <- function(
    X,
    y,
    g,
    sub_genes,
    parent_gene = NULL,
    depth = 0,
    max_depth = 5,
    min_node_size = 10
) {
  
  result <- list(
    gene_delta = numeric(0),
    tree_edges = data.frame(u = character(0), v = character(0))
  )
  
  if (depth >= max_depth) return(result)
  if (length(y) < 2 * min_node_size) return(result)
  if (length(unique(y)) == 1) return(result)
  
  if (is.null(parent_gene)) {
    candidate_genes <- sub_genes
  } else {
    nb <- neighbors(g, parent_gene)
    candidate_genes <- names(nb)
    candidate_genes <- intersect(candidate_genes, sub_genes)
  }
  
  candidate_genes <- intersect(candidate_genes, colnames(X))
  
  if (length(candidate_genes) == 0) return(result)
  
  sp <- best_split(
    X = X,
    y = y,
    candidate_genes = candidate_genes,
    min_node_size = min_node_size
  )
  
  if (is.null(sp$gene)) return(result)
  if (sp$delta <= 0) return(result)
  
  split_gene <- sp$gene
  
  result$gene_delta[split_gene] <- sp$delta
  
  if (!is.null(parent_gene)) {
    result$tree_edges <- rbind(
      result$tree_edges,
      data.frame(
        u = parent_gene,
        v = split_gene,
        stringsAsFactors = FALSE
      )
    )
  }
  
  left_id <- which(X[, split_gene] <= sp$cut)
  right_id <- which(X[, split_gene] > sp$cut)
  
  left_res <- build_graph_guided_tree(
    X = X[left_id, , drop = FALSE],
    y = y[left_id],
    g = g,
    sub_genes = sub_genes,
    parent_gene = split_gene,
    depth = depth + 1,
    max_depth = max_depth,
    min_node_size = min_node_size
  )
  
  right_res <- build_graph_guided_tree(
    X = X[right_id, , drop = FALSE],
    y = y[right_id],
    g = g,
    sub_genes = sub_genes,
    parent_gene = split_gene,
    depth = depth + 1,
    max_depth = max_depth,
    min_node_size = min_node_size
  )
  
  all_gene_delta <- c(result$gene_delta, left_res$gene_delta, right_res$gene_delta)
  gene_delta_sum <- tapply(all_gene_delta, names(all_gene_delta), sum)
  
  all_tree_edges <- rbind(
    result$tree_edges,
    left_res$tree_edges,
    right_res$tree_edges
  )
  
  return(list(
    gene_delta = gene_delta_sum,
    tree_edges = all_tree_edges
  ))
}
```

---


# Step 4: Build the Graph-Guided Random Forest

```r
Rtree <- 100
sub_size <- 50
max_depth <- 5
min_node_size <- 10

all_edges <- as_data_frame(g_prior, what = "edges")[, 1:2]
colnames(all_edges) <- c("u", "v")

edge_key <- function(u, v) {
  apply(
    cbind(u, v),
    1,
    function(z) paste(sort(z), collapse = "_")
  )
}

all_edges$key <- edge_key(all_edges$u, all_edges$v)

Iuv_vec <- setNames(rep(0, nrow(all_edges)), all_edges$key)

gene_importance_all <- setNames(rep(0, ncol(EXP_PMIEA)), colnames(EXP_PMIEA))

for (r in 1:Rtree) {
   
  seed_gene <- sample(V(g_prior)$name, 1)
  sub_genes <- get_bfs_subnetwork(g_prior, seed_gene, m = sub_size)
  
  if (length(sub_genes) < 2) next
  
  boot_id <- sample(seq_len(nrow(EXP_PMIEA)), replace = TRUE)
  
  X_boot <- EXP_PMIEA[boot_id, sub_genes, drop = FALSE]
  y_boot <- Y_PMIEA[boot_id]
  
  tree_res <- build_graph_guided_tree(
    X = X_boot,
    y = y_boot,
    g = g_prior,
    sub_genes = sub_genes,
    max_depth = max_depth,
    min_node_size = min_node_size
  )
  
  if (length(tree_res$gene_delta) == 0) next
  
  gene_importance_all[names(tree_res$gene_delta)] <-
    gene_importance_all[names(tree_res$gene_delta)] +
    tree_res$gene_delta
  
  if (nrow(tree_res$tree_edges) == 0) next
  
  tree_edges <- unique(tree_res$tree_edges)
  tree_edges$key <- edge_key(tree_edges$u, tree_edges$v)
  
  for (ee in seq_len(nrow(tree_edges))) {
    
    u <- tree_edges$u[ee]
    v <- tree_edges$v[ee]
    key <- tree_edges$key[ee]
    
    if (!key %in% names(Iuv_vec)) next
    
    delta_u <- ifelse(u %in% names(tree_res$gene_delta),
                      tree_res$gene_delta[u], 0)
    
    delta_v <- ifelse(v %in% names(tree_res$gene_delta),
                      tree_res$gene_delta[v], 0)
    
    Iuv_vec[key] <- Iuv_vec[key] + (delta_u + delta_v) / 2
  }
}
```

---


# Step 5: Compute Molecular Interplay Importance

```r
IuvSCORE <- all_edges
IuvSCORE$Iuv <- Iuv_vec[IuvSCORE$key]

IuvSCORE <- IuvSCORE[order(IuvSCORE$Iuv, decreasing = TRUE), ]

gene_importance_result <- data.frame(
  Gene = names(gene_importance_all),
  I_gene = as.numeric(gene_importance_all),
  stringsAsFactors = FALSE
)

gene_importance_result <- gene_importance_result[
  order(gene_importance_result$I_gene, decreasing = TRUE),
]

############################################################
## Normalize Iuv
############################################################

if(
  max(IuvSCORE$Iuv, na.rm = TRUE) ==
  min(IuvSCORE$Iuv, na.rm = TRUE)
){
  IuvSCORE$Iuv_norm <- 0
}else{
  IuvSCORE$Iuv_norm <-
    (
      IuvSCORE$Iuv -
      min(IuvSCORE$Iuv, na.rm = TRUE)
    ) /
    (
      max(IuvSCORE$Iuv, na.rm = TRUE) -
      min(IuvSCORE$Iuv, na.rm = TRUE)
    )
}
```

---


# Step 6: Compute the Topology-Aware Molecular Interplay Score

```r
############################################################
## Differential edge betweenness dB_uv
############################################################

tumor_edges <- as_data_frame(g_prior, what = "edges")[, 1:2]
colnames(tumor_edges) <- c("u", "v")
tumor_edges$key <- edge_key(tumor_edges$u, tumor_edges$v)

tumor_edges$B_tumor <- edge_betweenness(
  g_prior,
  directed = FALSE,
  weights = NA
)

EDGE_NORMAL <- NP_NM[, 1:2, drop = FALSE]
EDGE_NORMAL <- as.data.frame(EDGE_NORMAL)
colnames(EDGE_NORMAL) <- c("u", "v")

EDGE_NORMAL$u <- as.character(EDGE_NORMAL$u)
EDGE_NORMAL$v <- as.character(EDGE_NORMAL$v)

EDGE_NORMAL <- EDGE_NORMAL[
  EDGE_NORMAL$u %in% common_genes &
    EDGE_NORMAL$v %in% common_genes,
]

g_normal <- graph_from_data_frame(
  EDGE_NORMAL,
  directed = FALSE,
  vertices = common_genes
)

g_normal <- simplify(g_normal)

normal_edges <- as_data_frame(g_normal, what = "edges")[, 1:2]
colnames(normal_edges) <- c("u", "v")
normal_edges$key <- edge_key(normal_edges$u, normal_edges$v)

normal_edges$B_normal <- edge_betweenness(
  g_normal,
  directed = FALSE,
  weights = NA
)

edge_topology <- merge(
  tumor_edges[, c("key", "u", "v", "B_tumor")],
  normal_edges[, c("key", "B_normal")],
  by = "key",
  all = TRUE
)

edge_topology$B_tumor[is.na(edge_topology$B_tumor)] <- 0
edge_topology$B_normal[is.na(edge_topology$B_normal)] <- 0

missing_uv <- is.na(edge_topology$u) | is.na(edge_topology$v)

if (any(missing_uv)) {
  uv_mat <- do.call(
    rbind,
    strsplit(edge_topology$key[missing_uv], "_")
  )
  edge_topology$u[missing_uv] <- uv_mat[, 1]
  edge_topology$v[missing_uv] <- uv_mat[, 2]
}

edge_topology$dB_uv <- abs(
  edge_topology$B_tumor - edge_topology$B_normal
)

############################################################
## Differential hubness dH_uv
############################################################

deg_tumor <- degree(g_prior, mode = "all")
deg_normal <- degree(g_normal, mode = "all")

all_genes <- union(names(deg_tumor), names(deg_normal))

deg_tumor_all <- setNames(rep(0, length(all_genes)), all_genes)
deg_normal_all <- setNames(rep(0, length(all_genes)), all_genes)

deg_tumor_all[names(deg_tumor)] <- deg_tumor
deg_normal_all[names(deg_normal)] <- deg_normal

edge_topology$d_u_tumor <- deg_tumor_all[edge_topology$u]
edge_topology$d_v_tumor <- deg_tumor_all[edge_topology$v]

edge_topology$d_u_normal <- deg_normal_all[edge_topology$u]
edge_topology$d_v_normal <- deg_normal_all[edge_topology$v]

edge_topology$H_tumor <-
  edge_topology$d_u_tumor * edge_topology$d_v_tumor

edge_topology$H_normal <-
  edge_topology$d_u_normal * edge_topology$d_v_normal

edge_topology$dH_uv <- abs(
  edge_topology$H_tumor - edge_topology$H_normal
)

############################################################
## Merge differential topology measures with IuvSCORE
############################################################

IuvSCORE <- merge(
  IuvSCORE,
  edge_topology[, c(
    "key",
    "B_tumor", "B_normal", "dB_uv",
    "H_tumor", "H_normal", "dH_uv",
    "d_u_tumor", "d_v_tumor",
    "d_u_normal", "d_v_normal"
  )],
  by = "key",
  all.x = TRUE
)

IuvSCORE$dB_uv[is.na(IuvSCORE$dB_uv)] <- 0
IuvSCORE$dH_uv[is.na(IuvSCORE$dH_uv)] <- 0

############################################################
## Min-max normalization
############################################################

minmax_norm <- function(x) {
  x <- as.numeric(x)
  
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  
  if (xmax == xmin) {
    return(rep(0, length(x)))
  }
  
  return((x - xmin) / (xmax - xmin))
}

IuvSCORE$dB_uv_tilde <- minmax_norm(IuvSCORE$dB_uv)
IuvSCORE$dH_uv_tilde <- minmax_norm(IuvSCORE$dH_uv)

############################################################
## Topology-Aware Molecular Interplay Predictive Score
############################################################
IuvSCORE$W_uv <- 1 / (
  1 + exp(-(IuvSCORE$dB_uv_tilde +IuvSCORE$dH_uv_tilde+IuvSCORE$Iuv_norm))
)
```

---


# Step 7: Perform Pathway Enrichment Analysis

```r
# Number of permutation: 200
NOpm<-201 

SC_T<-matrix(numeric(1*(NOpm)),ncol=1)
for (pm in 1:(NOpm)){

if (pm==1){
Pathway_Genes<-Observed_Pathway_Genes
} else if (pm>1){
Pathway_Genes<- sample(unique(c(IuvSCORE[,"u"],IuvSCORE[,"v"])),length(Observed_Pathway_Genes))
}
EnSCORE<-matrix(numeric(length(IuvSCORE[,1])*6),ncol=6)
rownames(EnSCORE)<-IuvSCORE[,1]
colnames(EnSCORE)<-c("R","IN","NiN","Hit","Miss","SC")
EnSCORE[,1]<-IuvSCORE[,"W_uv"]
EnSCORE<-EnSCORE[order(EnSCORE[,1],decreasing=TRUE),]

EnSCORE0<-EnSCORE
EnSCORE0[is.na(match(read.csv(text=rownames(EnSCORE), sep="_", header=FALSE)[,1],Pathway_Genes))==FALSE,2]<-1
EnSCORE0[is.na(match(read.csv(text=rownames(EnSCORE), sep="_", header=FALSE)[,2],Pathway_Genes))==FALSE,3]<-1
EnSCORE[apply(EnSCORE0[,2:3],1,sum)==2,2]<-1
EnSCORE[EnSCORE[,2]==0,3]<-1

N<-nrow(EnSCORE)
Nh<-sum(EnSCORE[,2])
NR<-sum(EnSCORE[,2]) 
if (NR>0){EnSCORE[,4]<-cumsum(EnSCORE[,2])/NR}
if ((N-Nh)>0){EnSCORE[,5]<-cumsum(EnSCORE[,3])/(N-Nh)}
EnSCORE[,6]<-EnSCORE[,4]-EnSCORE[,5]

SC_T[pm ,1]<-max(EnSCORE[,6])
}
SC_T0<-SC_T
nmSC_T<-SC_T

nmSC_T[SC_T0[,1]>0,1]<-SC_T0[SC_T0[,1]>0,1]/mean(SC_T0[2:NOpm,][SC_T0[2:NOpm,1]>0])
nmSC_T[SC_T0[,1]<0,1]<-SC_T0[SC_T0[,1]<0,1]/mean(SC_T0[2:NOpm,][SC_T0[2:NOpm,1]<0])
Pvalue<-sum(nmSC_T[1,1]<=nmSC_T[-1,1][nmSC_T[-1,1]>0])/(NOpm-1)

print(Pvalue)
```

---


## Main Output Objects

| Object | Description |
|---|---|
| `gene_importance_result` | Phenotype-predictive importance of individual genes |
| `IuvSCORE` | Molecular interplay importance, differential topology measures, and topology-aware scores |
| `SC_T` | Observed and permutation-based enrichment scores |
| `nmSC_T` | Normalized enrichment scores |
| `Pvalue` | Empirical permutation p-value |

---

## Reproducibility

The analysis includes random seed-gene selection, bootstrap sampling, and gene-set permutation. Add the following line before the random forest analysis to obtain reproducible results:

```r
set.seed(1234)
```

---

## Citation

If you use this repository, please cite:

> Park H, Imoto S.  
> **PMIEA: Phenotype-Predictive Molecular Interplays Enrichment Analysis using Graph-Guided Random Forest.**  
> Manuscript in preparation.

---

## Contact

**Heewon Park**  
School of Mathematics, Statistics and Data Science  
Sungshin Women's University  
Email: heewonn.park@gmail.com
