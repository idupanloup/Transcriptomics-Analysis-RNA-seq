# --------------------------------
# Loading R packages

library(pheatmap)

library(clusterProfiler)
library(DESeq2)
library(org.Mm.eg.db)



# --------------------------------
# Data description and importation

data <- read.table("RNASEQ22_Day3_dataset_genecounts.txt", header = TRUE, row.names = 1)
data <- as.matrix(data)
meta <- read.table("RNASEQ22_Day3_dataset_infossamples.txt", header = TRUE)
rownames(meta) <- meta$labels



# --------------------------------
# Count normalization using DESeq2

# Match the metadata and counts data

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# Create DESEq2 object

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ group)

# Generate the normalized counts

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)



# --------------------------------
# QC for DE analysis using DESeq2

# Transform normalized counts using the rlog function

rld <- rlog(dds, blind=TRUE)

# Principal components analysis (PCA)

plotPCA(rld, intgroup="group")

# Hierarchical Clustering

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor)



# --------------------------------
# Differential expression analysis with DESeq2

# Run the differential expression analysis

dds <- DESeq(dds)

# Plot dispersion

plotDispEsts(dds)

# Build results table

res <- results(dds)
summary(res)

# plotMA

plotMA(res, ylim=c(-2,2))



# --------------------------------
# Functional analysis with clusterProfiler

# Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes <- as.character(rownames(res))

# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_genes <- as.character(rownames(signif_res))

# Run GO enrichment analysis 
ego <- enrichGO(gene = signif_genes, 
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

# Dotplot 
dotplot(ego, showCategory=50)

# Enrichmap
emapplot(ego, showCategory = 50)

# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange

# Category Netplot
cnetplot(ego, 
         categorySize ="pvalue", 
         showCategory = 5, 
         foldChange = signif_res_lFC, 
         vertex.label.font=6)

