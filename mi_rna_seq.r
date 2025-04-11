
library(DESeq2)
library(ggplot2)
library(dplyr)
library(randomForest)
library(clusterProfiler)
library(org.Hs.eg.db)

# Simulated miRNA count matrix
table <- matrix(sample(20:1000, 40, replace=TRUE), nrow=10)
rownames(table) <- paste0("miR-", 1:10)
colnames(table) <- c("tumor1", "tumor2", "normal1", "normal2")

meta <- data.frame(
  condition = factor(c("tumor", "tumor", "normal", "normal")),
  row.names = colnames(table)
)

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = table, colData = meta, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Volcano plot
res$miRNA <- rownames(res)
res_df <- as.data.frame(res)
res_df <- res_df %>% mutate(sig = padj < 0.1)

volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue), color=sig)) +
  geom_point() + theme_minimal() + ggtitle("miRNA Differential Expression")
print(volcano)

# Barplot of top miRNAs
top_miRNAs <- res_df %>% arrange(padj) %>% head(5)
ggplot(top_miRNAs, aes(x=reorder(miRNA, -baseMean), y=baseMean)) +
  geom_bar(stat='identity', fill='steelblue') +
  theme_minimal() +
  ggtitle("Top Expressed miRNAs") +
  ylab("Base Mean")

# Save top miRNAs to CSV

top_hits <- res_df %>% arrange(padj) %>% dplyr::select(miRNA, padj, log2FoldChange, baseMean) %>% head(10)

write.csv(top_hits, "top_miRNA_hits.csv", row.names = FALSE)



# Pathway annotation for top DE miRNAs (convert to ENTREZ if known)
# This part is illustrative – real use would need actual miRNA → gene target mapping
# Here we simulate gene target IDs
example_genes <- sample(keys(org.Hs.eg.db, keytype="ENTREZID"), 50)
pathway_enrich <- enrichKEGG(gene = example_genes, organism = "hsa")

# Save pathway results
write.csv(as.data.frame(pathway_enrich), "mirna_pathway_enrichment.csv", row.names = FALSE)
