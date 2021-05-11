# load packages
setwd("DESeq2//")
library(apeglm)
library(DESeq2)
library(pheatmap)
library(plotly)
library(data.table)

# import data
countdata <- read.table("FeatureCounts.csv", header = T)
gene_id <- read.csv("gene_id_L1HS.txt")
gene_name <- read.csv("gene_name_L1HS.txt")
gene_type <- read.csv("gene_type_L1HS.txt")


## Data preparation
#---------------------------------------------------------------------
# how many lines 
nrow(gene_id)
nrow(gene_name)
nrow(gene_type)


# concatinate gene_id, gene_name, gene_type and make unique
gene_info <- data.frame (gene_id, gene_name, gene_type)
gene_info <- unique(gene_info)
nrow(gene_info )
test <- unique(countdata$Geneid)

# Delete and rename columns
countdata <- data.frame(countdata$Geneid, countdata$...Results.Bulk.RNA.2_Mapping.Bulk.RNA_S1_B1_Aligned.sortedByCoord.out.bam, countdata$...Results.Bulk.RNA.2_Mapping.Bulk.RNA_S1_B2_Aligned.sortedByCoord.out.bam, countdata$...Results.Bulk.RNA.2_Mapping.Bulk.RNA_S2_B1_Aligned.sortedByCoord.out.bam, countdata$...Results.Bulk.RNA.2_Mapping.Bulk.RNA_S2_B2_Aligned.sortedByCoord.out.bam, countdata$...Results.Bulk.RNA.3_Chimp_Mapping.S1_B1_Aligned.sortedByCoord.out.bam, countdata$...Results.Bulk.RNA.3_Chimp_Mapping.S1_B2_Aligned.sortedByCoord.out.bam, countdata$...Results.Bulk.RNA.3_Chimp_Mapping.S2_B1_Aligned.sortedByCoord.out.bam, countdata$...Results.Bulk.RNA.3_Chimp_Mapping.S2_B2_Aligned.sortedByCoord.out.bam)
Newnames <- c('GeneID','H1_B1', 'H1_B2', 'H2_B1', 'H2_B2', 'Ch1_B1', 'Ch1_B2', 'Ch2_B1', 'Ch2_B2')
names(countdata) <- Newnames

# Rename gene info and make it unique
info.names <- c('GeneID', 'GeneName', 'GeneType')
names(gene_info) <- info.names

#sort both files before merge together
gene_info <- gene_info[order(gene_info$GeneID),]
countdata <- countdata[order(countdata$GeneID),]
# check if identical gene_id
all(gene_info$gene_id == countdata$gene_id)

# merge countdata and gene_info
countdata_gene <- merge(gene_info,countdata, by="GeneID")

# gene_name as rownames
rownames(countdata_gene) <- make.unique(as.character(countdata_gene$GeneName))

# delete 3 columns
countdata_gene <- countdata_gene[-1]
countdata_gene <- countdata_gene[-1]
countdata_gene <- countdata_gene[-1]
countdata <- countdata_gene


############################################################################
# *Metadata
(condition <- factor(c("human", "human", "human", "human", "chimp", "chimp", "chimp", "chimp")))
(coldata <- data.frame(row.names=colnames(countdata),sample=colnames(countdata), condition))
#------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
## DESeq2 analysis

# Make sure the counts are numeric
sapply(countdata, class)
countdata <- transform(countdata, H1_B1=as.numeric(H1_B1), H1_B2=as.numeric(H1_B2), H2_B1=as.numeric(H2_B1), H2_B2=as.numeric(H2_B2), Ch1_B1=as.numeric(Ch1_B1), Ch1_B2=as.numeric(Ch1_B2), Ch2_B1=as.numeric(Ch2_B1), Ch2_B2=as.numeric(Ch2_B2))
sapply(countdata, mode)

# DESeq object
dds <- DESeqDataSetFromMatrix (countData = countdata, colData = coldata, design = ~condition)

# Filter out lowly expressed counts
dds <- dds[rowSums(counts(dds))>0.5,]

# *Set factor and run DESeq [Genes that are upreg in human when compared to chimp]
dds$condition <- factor(dds$condition, levels = c('chimp', 'human'))
dds <- DESeq(dds)

# DESeq results object
res <- results(dds)
summary(res)
res

# Regularized log transformation for visualization
# Normalization
norm <- counts(dds, normalized=TRUE)
norm_dds <- norm
#rlogTransformation transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size
rld <- rlogTransformation(dds)

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_human_vs_chimp", type="apeglm")

# Adjusted p-value, up and down-regulated genes
table(res$padj < 0.05)
table(res$padj < 0.05 & res$log2FoldChange > 2)
table(res$padj < 0.05 & res$log2FoldChange < -2)


# Print to table expressed, up and downregulated genes
up_regulated_L1HS <- subset(res, res$padj < 0.05 & res$log2FoldChange > 2)
down_regulated_L1HS <- subset(res, res$padj < 0.05 & res$log2FoldChange < -2)
all_L1HS <- subset(res)

write.table(up_regulated_L1HS, file="L1HS_HumvsChim\\up_regulated_L1HS.txt", sep = "\t", quote = FALSE)
write.table(down_regulated_L1HS, file= "L1HS_HumvsChim\\down_regulated_L1HS.txt", sep = "\t", quote = FALSE)
write.table(all_L1HS, file="L1HS_HumvsChim\\all_regulated_L1HS.txt", sep = "\t", quote = FALSE)

## Visualization
# Principal Component Analysis
png("L1HS_HumvsChim\\PCA_L1HS.png")
plotPCA(rld) + theme_classic()
dev.off()
# MA plot
plotMA(res, alpha=0.05, ylim=c(-6,6))

# MA plot with shrunken log2 fold changes.
png("L1HS_HumvsChim\\MA_L1HS.png")
plotMA(resLFC, ylim=c(-4,4))
dev.off()
# Mean plot
norm_dds <- as.data.frame(norm_dds)
norm_dds_test <- norm_dds
# Mean of human and chimp
human_mean <- rowMeans(norm_dds_test[,as.character(subset(coldata, coldata$condition == "human")$sample)])
chimp_mean <- rowMeans(norm_dds_test[,as.character(subset(coldata, coldata$condition == "chimp")$sample)])

# Genes found in humans and chimp
human_genes_expressed <- norm_dds_test[,subset(coldata, coldata$condition == "Human")$sample]
human_genes_expressed <- human_genes_expressed[rowSums(human_genes_expressed)>0,]
chimp_genes_expressed <- norm_dds_test[,subset(coldata, coldata$condition == "Chimp")$sample]
chimp_genes_expressed <- chimp_genes_expressed[rowSums(chimp_genes_expressed)>0,]

# Mark up-regulated genes red and down-regulated genes blue
norm_dds_test$cex <- rep(1, nrow(norm_dds_test))
norm_dds_test[rownames(up_regulated_L1HS),] <- 2
norm_dds_test[rownames(down_regulated_L1HS),] <- 3
norm_dds_test$colour <- "black"
norm_dds_test$colour[norm_dds_test$cex==2] <- "firebrick"
norm_dds_test$colour[norm_dds_test$cex==3] <- "steelblue"
# Plot with scatter plot
png("L1HS_HumvsChim\\Thesis\\MeanScatter_L1HS.png")
plot(log2(chimp_mean + 0.5),
     log2(human_mean + 0.5),
     col=norm_dds_test$colour,
     pch=16,
     xlab = "log2(mean Chimp)",
     ylab = "log2(mean Human)")
legend("bottomright",c("upregulated","not significant", "downregulated"),cex= .8, fill=c("firebrick", "black", "steelblue"))
dev.off()

# Heatmaps
# select condition to be compared in heatmaps
coldata <- coldata[order(coldata$condition),]
# up and down-regulated genes
signdiff_up <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange > 2)),]
signdiff_low <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange < -2)),]

# Heatmap of top variance of genes
topVarGenes <- head(order(-rowVars(assay(rld))),20)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("condition"),drop=FALSE])
png("L1HS_HumvsChim\\topHeat_L1HS.png")
pheatmap(mat, annotation_col = df)
dev.off()
# Heatmap of upregulation
png("L1HS_HumvsChim\\upHeat_L1HS.png")
pheatmap(log2(signdiff_up[,as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = T, annotation_col = df, fontsize = 8)
dev.off()
# Heatmap of downregulation
png("L1HS_HumvsChim\\downHeat_L1HS.png")
pheatmap(log2(signdiff_low[,as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = T, annotation_col = df, fontsize = 8)
dev.off()