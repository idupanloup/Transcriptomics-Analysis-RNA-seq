# --------------------------------
# Loading R packages

library(devtools)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)

library(DESeq2)
library(edgeR)
library(mixOmics)



# --------------------------------
# Data description and importation

raw_counts <- read.table("RNASEQ22_Day2_dataset_genecounts.txt", header = TRUE, 
                         row.names = 1)
raw_counts <- as.matrix(raw_counts)
gene_lengths <- scan("RNASEQ22_Day2_dataset_genelength.txt")
design <- read.table("RNASEQ22_Day2_dataset_infossamples.txt", header = TRUE)



# --------------------------------
# Basic exploratory analysis of raw counts

raw_counts_wn <- raw_counts[rowSums(raw_counts) > 0, ]
dim(raw_counts_wn)

log_counts <- log2(raw_counts_wn + 1)
head(log_counts)

df_raw <- melt(log_counts, id = rownames(raw_counts_wn))
names(df_raw)[1:2]<- c("id", "sample")
df_raw$method <- rep("Raw counts", nrow(df_raw))  
head(df_raw)

# Count distribution

df <- data.frame(rcounts = raw_counts_wn [ ,1], lcounts = log_counts[ ,1])

p1 <- ggplot(data=df, aes(x = rcounts, y = ..density..))
p1 <- p1 + geom_histogram(fill = "lightblue")
p1 <- p1 + theme_bw()
p1 <- p1 + ggtitle(paste0("count distribution '", design$labels[1], "'"))
p1 <- p1 + xlab("counts")

p2 <- ggplot(data=df, aes(x = lcounts, y = ..density..))
p2 <- p2 + geom_histogram(fill = "lightblue")
p2 <- p2 + theme_bw()
p2 <- p2 + ggtitle(paste0("count distribution - '", design$labels[1], "'"))
p2 <- p2 + xlab(expression(log[2](counts + 1)))

grid.arrange(p1, p2, ncol = 2)

# Relation between mean and variance

df <- data.frame(mean = apply(raw_counts_wn[ ,design$group == "newborns"], 1, mean),
                 var = apply(raw_counts_wn[ ,design$group == "newborns"], 1, var))
df <- df[df$mean <= 5000, ]
p <- ggplot(data=df, aes(x = mean, y = var))
p <- p + geom_point(colour = "orange")
p <- p + theme_bw()
p <- p + geom_abline(aes(intercept=0, slope=1))
p <- p + ggtitle("Variance versus mean in counts") + ylab("variance")
print(p)



# --------------------------------
# Normalization

# DESeq

groups <- factor(design$group)  
# DESeqDataSetFromMatrix(countData,colData,design)
# countData: a matrix of non-negative integers
# colData: a data.frame with at least a single column. Rows of colData correspond to columns of countData
# design: a formula or matrix: the formula expresses how the counts for each gene depend on the variables in colData.
dds <- DESeqDataSetFromMatrix(raw_counts_wn, DataFrame(groups), ~ groups)
dds
head(counts(dds))

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

deseq_normcount <- counts(dds, normalized = TRUE)

pseudo_deseq <- log2(deseq_normcount + 1)
df_deseq <- melt(pseudo_deseq, id = rownames(raw_counts_wn))
names(df_deseq)[1:2]<- c("id", "sample")
df_deseq$method <- rep("DESeq", nrow(df_raw))

# edgeR

dge2 <- DGEList(raw_counts_wn)
dge2

# calculate scaling factors to convert raw library sizes into effective library sizes
dge2 <- calcNormFactors(dge2, method = "TMM")

pseudo_TMM <- log2(cpm(dge2) + 1)

df_TMM <- melt(pseudo_TMM, id = rownames(raw_counts_wn))
names(df_TMM)[1:2] <- c ("id", "sample")
df_TMM$method <- rep("TMM", nrow(df_TMM))

# RPKM

gene_lengths_wn <- gene_lengths[rowSums(raw_counts) > 0]
pseudo_RPKM <- log2(rpkm(dge2, gene.length = gene_lengths_wn) + 1)

df_RPKM <- melt(pseudo_RPKM, id = rownames(raw_counts_wn))
names(df_RPKM)[1:2] <- c ("id", "sample")
df_RPKM$method <- rep("RPKM", nrow(df_RPKM))

# Comparison

df_allnorm <- rbind(df_raw, df_deseq, df_TMM, df_RPKM)
df_allnorm$method <- factor(df_allnorm$method,
                            levels = c("Raw counts", "DESeq", "TMM", "RPKM"))

p <- ggplot(data=df_allnorm, aes(x=sample, y=value, fill=method))
p <- p + geom_boxplot()  
p <- p + theme_bw()
p <- p + ggtitle("Boxplots of normalized pseudo counts\n
for all samples by normalization methods")
p <- p + facet_grid(. ~ method) 
p <- p + ylab(expression(log[2] ~ (normalized ~ count + 1))) + xlab("")
p <- p + theme(title = element_text(size=10), axis.text.x = element_blank(), 
               axis.ticks.x = element_blank())
print(p)

