---
title: Transcriptomics Analysis RNA-seq
summary: course website
author: Isabelle Dupanloup
date: 2022-01-06
some_url: https://idupanloup.github.io/Transcriptomics-Analysis-RNA-seq/
---

---------------------
Knowledge / competencies

- Participants should already have a basic knowledge of Next Generation Sequencing (NGS) techniques; this course will discuss only the data analysis steps and not the data generation.
- A basic knowledge in statistics is required. Participants should know about p-values, student T-test, multiple testing correction and classification, PCA.
- A basic knowledge of R is also required. Participants should know how to read files, run PCA, do classification, visualise heatmaps using R command lines.

---------------------
Requirements

- Hardware (64-bit computer with 4 GB of RAM (8 GB preferred))

- FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- QualiMap (http://qualimap.bioinfo.cipf.es/)

- R (https://www.r-project.org, version > 4.0)
- latest version of R Studio

---------------------
Installation of R packages for Practical 1

- if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

- BiocManager::install("NOISeq")
- BiocManager::install("Repitools")
- BiocManager::install("Rsamtools")
- BiocManager::install("Rsubread")
- BiocManager::install("rtracklayer")

---------------------
Installation of R packages for Practical 2

- install.packages("devtools")
- install.packages("ggplot2")
- install.packages("gridExtra")
- install.packages("RColorBrewer")
- install.packages("reshape2")

- if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

- BiocManager::install("DESeq2")
- BiocManager::install("edgeR")
- BiocManager::install("mixOmics")

---------------------
Installation of R packages for Practical 3

- install.packages("pheatmap")

- if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

- BiocManager::install("clusterProfiler")
- BiocManager::install("org.Mm.eg.db")
