BiocManager::install("DESeq2")
library(tidyverse, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(biomaRt, quietly = TRUE)
sampleGroup <- read.table("sampleGroup.csv", header=T)
head(sampleGroup)
sampleGroup$Group <- factor(sampleGroup$Group, levels = c("Polytoe", "normal"))
print(sampleGroup)
readCount <-  read.table("gene_count_matrix.csv",sep=",",header=T)
rownames(readCount)=readCount[,1]
readCount <- readCount[,-1]
head(readCount)
readCount <- as.data.frame(readCount)
readCount[] <- lapply(readCount, as.integer)
print(head(readCount, n = 3))
dim(readCount)
dim(sampleGroup)
colnames(readCount)
rownames(sampleGroup)
dds <- DESeqDataSetFromMatrix(countData = readCount, colData = sampleGroup, design = ~Group)
dds
dds_mean <- DESeq(dds, fitType = "mean", minReplicatesForReplace = 7)
plotDispEsts(dds_mean)
dds_param <- DESeq(dds, fitType = "parametric", minReplicatesForReplace = 7)
plotDispEsts(dds_param)
dds1 <- DESeq(dds, fitType = 'parametric', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c("Group", "Polytoe", "normal"))
res
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'Polytoe_normal.DESeq2.csv', col.names = NA, sep = '\t', quote = FALSE)
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'none'

write.table(res1, file = 'Polytoe_normal.DESeq2.all.csv', sep = '\t', col.names = NA, quote = FALSE)
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'Polytoe_normal.DESeq2.select.csv', sep = '\t', col.names = NA, quote = FALSE)

res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'Polytoe_normal.DESeq2.up.csv', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'Polytoe_normal.DESeq2.down.csv', sep = '\t', col.names = NA, quote = FALSE)
