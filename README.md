# eQTL-Based Gene Regulatory Network Inference Pipeline

This repository contains a complete workflow for processing gene expression data, mapping expression quantitative trait loci (eQTLs), prioritizing candidate regulatory variants, and inferring gene regulatory networks (GRNs) using eQTLnet.

## Workflow Overview

The pipeline consists of the following steps:

### 1. Expression Matrix Preparation
**Script:** `01_prepare_expression_matrices`

- Import raw expression data.
- Format expression matrices for downstream analyses.
- Harmonize sample and gene identifiers.

### 2. Low Expression Filtering
**Script:** `02_filter_low_expression_genes`

- Remove genes with insufficient expression levels.
- Reduce noise and improve statistical power.

### 3. Removal of Replicate-Specific Zero Expression
**Script:** `03_remove_replicate_specific_zeros`

- Identify genes with expression detected in only a subset of biological replicates.
- Remove genes with unreliable expression profiles.

### 4. Estimation of Expression BLUPs
**Script:** `04_estimate_expression_blups`

- Estimate Best Linear Unbiased Predictors (BLUPs) for gene expression.

### 5. eQTL Mapping
**Script:** `05_eQTL_mapping_MatrixEQTL`

- Perform genome-wide eQTL mapping using MatrixEQTL.

### 6. SNP Prioritization
**Script:** `06_eQTL_SNP_Prioritization`

- Prioritize significant eQTL SNPs.
- Select representative markers for downstream analyses.

### 7. eQTL Gene Set Construction
**Script:** `07_eQTL_Gene_Set`

- Generate gene sets associated with significant eQTL-genes

### 8. Identification of Genes in eQTL Intervals
**Script:** `08_Identify_Genes_in_eQTL_Intervals`

- Identify genes physically located within eQTL confidence intervals.

### 9. Preparation of eQTLnet Inputs
**Script:** `09_Prepare_eQTLnet_inputs`

- Format genotype, expression, and candidate regulator data.
- Generate input files required by eQTLnet.

### 10. eQTLnet Regression Analysis
**Script:** `10_eQTLnet_regression`

- Run eQTLnet regression analysis.
- Infer putative regulatory relationships.

### 11. Construction of Network Matrices
**Script:** `11_Build_Beta_Network_Matrices`

- Convert eQTLnet outputs beta coefficient matrices.
- Generate network representations suitable for visualization and analysis.

### 12. Benchmarking of Inferred Networks
**Script:** `12_Benchmark_GRN_Against_GS`

- Compare inferred gene regulatory networks against reference gene sets.
- Evaluate network performance and biological relevance.

---

## Dependencies

Main R packages used throughout the pipeline include:

- MatrixEQTL
- data.table
- dplyr
- tidyr
- readr
- GenomicRanges
- eQTLnet
- igraph

Additional package requirements are specified within individual scripts.

---

## Input Data

The pipeline requires:

- Gene expression matrix
- Genotype matrix
- Marker position information
- Gene annotation files
- Experimental metadata

---
