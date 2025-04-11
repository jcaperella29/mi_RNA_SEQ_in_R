Description
Runs a full miRNA-seq pipeline: differential expression, and pathway enrichment analysis.

Features

DESeq2-based differential expression analysis

Volcano and bar plots of top miRNAs


KEGG pathway enrichment via clusterProfiler

Output CSV files:

top_miRNA_hits.csv

mirna_pathway_enrichment.csv

Requirements

R packages: DESeq2, ggplot2, dplyr, randomForest, clusterProfiler, org.Hs.eg.db

Usage
source("mirna_seq_script.R")


