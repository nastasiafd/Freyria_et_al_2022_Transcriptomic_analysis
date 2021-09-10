######### Scripts for "Salinity tolerance and stress responses in an Arctic Pelagophyte revealed by comparative transcriptomic and gene expression analysis" 
######### by Nastasia J. Freyria, Alan Kuo, Mansi Chovatia, Jenifer Johnson, Anna Lipzen, Kerrie W. Barry, Igor V. Grigoriev and Connie Lovejoy.

# Computer code is following methods from packages DESeq2 by Michael I. Love, Simon Anders, and Wolfgang Huber
# and DESeq by Simon Anders

# Libraries
library(DESeq2) # v.1.26.0
library(DESeq) # v.1.38.0 
library(vegan) # v.2.5-7
library(gplots) # v.3.1.1
library(ggplot2) # v.3.3.5
library(pheatmap) # v.1.0.12
library(UpSetR) # v.1.4.0


####### Create matrix DESEQ2 ####### 

# Read files
Transcrip_all_counts = read.table("all_transcriptomes_featurecounts.txt", row.names = 1, header = T)

# Define replicates for both conditions to be compared
condition = factor(c(rep("45psu", 3), rep("35psu", 3), rep("25psu", 3), rep("16psu", 3), rep("8psu", 2)))

# Create a countDataSet
dds = DESeqDataSetFromMatrix(Transcrip_all_counts, DataFrame(condition), design = ~condition)

# Estimate normalization factors
dds = estimateSizeFactors(dds)
sizeFactors(dds)

# Estimate dispersion
dds = estimateDispersions(dds)

# Differential analysis
dds  = nbinomWaldTest(dds)
Results = results(dds)

# Result writing
write.table(Results, file = "DESeq2_results_all_transcriptome.txt", sep = "\t", quote = F)



####### Matrix analysis ####### 

# Import DESeq2 results
DESeq2Results = as.matrix(read.table("DESeq2_results_all_transcriptome.txt", header = T, row.names = 1))

# Vector with conditions
conditions_2 = factor(c(rep("D0", 2), rep("D1", 3), rep("D2", 3), rep("D3", 3), rep("D4", 3), rep("45psu", 3), rep("35psu", 3), rep("25psu", 3), rep("16psu", 3), rep("8psu", 2)))

# Create a countDataSet
dds = DESeqDataSetFromMatrix(Transcrip_all_counts, DataFrame(conditions_2), design = ~conditions_2)

# Graphics to observe before normalization
barplot(colSums(counts(dds)), col = "light blue", ylab = "Total number of read counts")

# Estimate normalization factors
dds = estimateSizeFactors(dds)
sizeFactors(dds)

# Estimate data reproducibility
dds = estimateDispersions(dds)
plotDispEsts(dds)

# p-values distribution
hist(DESeq2Res[,"padj"], nclass = 100, col = "light blue", xlab = "p-values calculated", main = "Distribution of adjusted p-values")

# Table with count values (after normalization)
dataRNAseq2 = counts(dds, normalized = TRUE)

# Graphics to observe after normalization
barplot(colSums(dataRNAseq2), col = "light blue", ylab = "Total number of read counts")

# Mean values for conditions replicates
meanT1C45counts = apply(dataRNAseq2[,1:2], 1, mean) # 45psu
meanT2C5counts = apply(dataRNAseq2[,3:5], 1, mean) # 45psu
meanT3C5counts = apply(dataRNAseq2[,6:8], 1, mean) # 45psu
meanT4C5counts = apply(dataRNAseq2[,9:11], 1, mean) # 45psu
meaNT5C45counts = apply(dataRNAseq2[,12:14], 1, mean) # 45psu
mean45counts = apply(dataRNAseq2[,15:17], 1, mean) # 45psu
mean35counts = apply(dataRNAseq2[,18:20], 1, mean) # 35psu
mean25counts = apply(dataRNAseq2[,21:23], 1, mean) # 25psu
mean16counts = apply(dataRNAseq2[,24:26], 1, mean) # 16psu
mean8counts = apply(dataRNAseq2[,27:28], 1, mean) # 8psu

# Log Fold Change 
logFC_T1C45_T2C45 = log2((meanT1C45counts + 1)/(meanT2C5counts + 1))
logFC_T2C45_T3C45 = log2((meanT2C5counts + 1)/(meanT3C5counts + 1))
logFC_T3C45_T4C45 = log2((meanT3C5counts + 1)/(meanT4C5counts + 1))
logFC_T4C45_T5C45 = log2((meanT4C5counts + 1)/(meaNT5C45counts + 1))
logFC_45_35 = log2((mean45counts + 1)/(mean35counts + 1))
logFC_45_25 = log2((mean45counts + 1)/(mean25counts + 1))
logFC_45_16 = log2((mean45counts + 1)/(mean16counts + 1))
logFC_45_8 = log2((mean45counts + 1)/(mean8counts + 1))
logFC_35_25 = log2((mean35counts + 1)/(mean25counts + 1))
logFC_35_16 = log2((mean35counts + 1)/(mean16counts + 1))
logFC_35_8 = log2((mean35counts + 1)/(mean8counts + 1))
logFC_25_16 = log2((mean25counts + 1)/(mean16counts + 1))
logFC_25_8 = log2((mean25counts + 1)/(mean8counts + 1))
logFC_16_8 = log2((mean16counts + 1)/(mean8counts + 1))

## Normalized data
# 45psu/35psu
hist(logFC_45_35, nclass = 100, main = "logFC (Pel_T_45psu/35psu) distribution \n(norm. data)", xlab = "log(45psu/35psu) value")
abline(v = 0, col = "forestgreen")
abline(v = 2, col = "deepskyblue", lty = "dashed")
abline(v = -2, col = "firebrick", lty = "dashed")

# Up- and down-regulated genes
upGenes_45_35 = names(logFC_45_35[logFC_45_35 > 2]) 
downGenes_45_35 = names(logFC_45_35[logFC_45_35 < -2])

# Gene expression level
exprLevel = apply(dataRNAseq2, 1, mean)

# logFC2 versus the level of gene expression --> Figure 3a and 4a 
plot(log(exprLevel), logFC_45_35, pch = 20, xlab = "", ylab = "", col = "gray", ylim = c(-10, 6), yaxt = 'n', xlim = c(-5, 15), xaxt = 'n', cex=1) +
  axis(2, seq(-10,6, by=2), col.axis="black", col = "black", las = 2, cex.axis=1.2) +
  axis(1, seq(0,10, by=5), col.axis="black", col = "black", las = 1, cex.axis=1.2) +
  mtext("Log gene expression level", cex=1.2, side=1, line=2.5, col="black") + 
  mtext("Log2 fold change", cex=1.2, side=2, line=2.5, col="black") +
  abline(h = 2, col = "black", lty = "dashed") + 
  abline(h = -2, col = "black", lty = "dashed") +
  points(log(exprLevel[upGenes_45_35]), logFC_45_35[upGenes_45_35], pch = 20, col = "firebrick2", cex=1.2) +
  points(log(exprLevel[downGenes_45_35]), logFC_45_35[downGenes_45_35], pch = 20, col = "darkorange1", cex=1.2)

# Do the same for all comparisons ...










####### Canonical Correspondence Analysis (CCA) --> Figure 1a ####### 

# Read file reads per gene 
T_read = read.csv("T_all_reads_per_genes.csv", row.names = 1, check.names = F)

# Filter by deleting zero values
T_reads = subset(T_read, rowSums(T_read)!=0)

# Transpose 
T_reads_t <- t(T_reads)

# Read file with metadata or environmental variables
env = read.csv("variables_reads.csv", row.names = 1, check.names = F)

# Verify if row names of both files are the same
env$Description == rownames(T_reads_t)

# If not the same order of labels, ordered it
env<-env[rownames(T_reads_t),]

# Verify again 
env$Description == rownames(T_reads_t)

# Calculate CCA  
T_reads.cca = cca(T_reads_t~., ENV, dist="bray", direction = 'forward', permutations = 999)

# Calculate adjusted R2 value
adjR2.ccaP <- RsquareAdj (T_reads.cca)$adj.r.squared 

# Incorpore adjustedR2 into the CCA calculation
T_reads.cca = cca(T_reads_t~., ENV, dist="bray", direction = 'forward', permutations = 999, R2scope = adjR2.ccaP)

# ANOVA test
anova(T_reads.cca)

# Calculate environmental variables significance
(fit.cca = envfit(T_reads.cca, ENV, perm=999, display="lc", scaling="sites", na.rm = TRUE))

# Plot CCA
plot(T_reads.cca, axes = T, ylab="", xlab="") +
  legend("bottomright",legend=c("R2adj = 0.716"), bty="n", cex=0.9) +
  plot(fit.cca, p.max=0.01, col="black", cex=0.9)



####### Principal Component Analysis (PCA) --> Figure 1b ####### 

# Define colors for conditions
cols <- c("red", "orange", "yellow", "navyblue", "turquoise")

# Calculated Regularized log transformation (rld)
rld <- rlog(dds, blind=TRUE)

# Plot PCA
plotPCA(rld, intgroup = c( "condition_2", "condition_2"))



####### Hierarchical clustering with all conditions --> Figure 1c ####### 

# Read transcriptome matrix calcultaed with featurecounts
tree = Transcrip_all_counts

# Define replicates for both conditions to be compared
conds = c(rep("D0", 2), rep("D1", 3), rep("D2", 3), rep("D3", 3), rep("D4", 3), rep("45psu", 3), rep("35psu", 3), rep("25psu", 3), rep("16psu", 3), rep("8psu", 2))

# Transform to data frame
expDesign = data.frame(row.names = colnames(tree), condition = conds, libType = rep("DESeq2_results_all_transcriptome", length(conds)))

# Create a CountDataSet object
cds = newCountDataSet(tree, expDesign$condition)

# Estimate size factor
cds = estimateSizeFactors(cds)
sizeFactors(cds)

# Estimate dispersion
cds = estimateDispersions(cds, method="blind", sharingMode="fit-only")
str(fitInfo(cds))

# Gene count in log
cnt = log(1 + counts(cds, normalized = T))

# Calculated distance
dst = as.dist(0.5 - 0.5 * cor(cnt))

#The complete linkage method finds similar clusters. 
plot(hclust(dst, method = "complete"))



####### Heatmap for Top25 genes --> Figure 2 ####### 

# Read file: matrix from FeatureCounts
Transcrip_all_counts = read.table("all_transcriptomes_featurecounts.txt", row.names = 1, header = T)

# Vector with conditions
condition_2 = as.factor(c(rep("D0", 2), rep("D1", 3), rep("D2", 3), rep("D3", 3), rep("D4", 3), rep("45psu", 3), rep("35psu", 3), rep("25psu", 3), rep("16psu", 3), rep("8psu", 2)))

# Create a countDataSet
dds = DESeqDataSetFromMatrix(Transcrip_all_counts, DataFrame(condition_2), design = ~condition_2)

# Estimate size factor
dds <- estimateSizeFactors(dds)

# Estimate dispersion
dds = estimateDispersions(dds)

# Normalization
dataRNAseq2 = counts(dds, normalized = TRUE)

# Expression gene level
exprLevel = apply(dataRNAseq2, 1, mean)

# highest variance between genes with vst (Variance stabilizing transformation)
vsd <- vst(dds, blind=FALSE)

# Selection of the top 25 genes
topVarGenes <- head(order(-rowVars(assay(vsd))),50)

# Create matrix
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat)

# Transform into a data frame
df <- as.data.frame(colData(vsd))

# Define conditions
ann_color = list(condition_2 = c(D0="darkred", D1="darkred", D2="darkred", D3="darkred", D4="darkred",
                                 `45psu`="firebrick2",`35psu`="darkorange", `25psu`="gold", `16psu`="dodgerblue2", `8psu`="turquoise2"))

# Pheatmap plot
pheatmap(mat, annotation_col=df, annotation_colors = ann_color,cutree_rows = 2, cutree_cols = 2)



####### Volcano plot --> Figure 3b and 4b #######

# Read file of two condition comparison
DESeq2_results_45_vs_8_logFC_45_vs_8 = read.table("DESeq2_results_45_vs_8_logFC_45_vs_8.txt", h=T, row.names = 1)

# Transform it into data frame
df = as.data.frame(DESeq2_results_45_vs_8_logFC_45_vs_8)

# Choose a threshold for p-value adjusted to highlight genes
up_adj = subset(df, abs(df$padj<0.05) & df$logFC45_8 > 2)
down_adj = subset(df, abs(df$padj<0.05) & df$logFC45_8 < -2)

# Choose a threshold for up- and down-regulated genes 
up_genes <- df$logFC35_25 > 2 & df$padj < 0.05
down_genes <- df$logFC35_25 < -2 & df$padj < 0.05

# Choose a threshold to separate up and down regulated genes on the plot
alpha <- 0.05 

# Volcano plot
plot(df$logFC45_8, -log10(df$padj), pch = 20, xlab = bquote(~Log[2]~fold~change), 
     ylab = bquote(~-log[10]~p~adj),
     main = "", col="gray", cex.axis=1.5,
     ylim = c(-2, 100), xlim = c(-10, 8),
     xaxt = 'n', yaxt = 'n') +
  abline(v = 2, col = "black", lty = "dashed") +
  abline(v = -2, col = "black", lty = "dashed") +
  abline(h=-log10(alpha), col="black", lty="dashed") +
  points(up_adj$logFC45_8, -log10(up_adj$padj), pch = 20, col = "darkorange1") +
  points(down_adj$logFC45_8, -log10(down_adj$padj), pch = 20, col = "goldenrod1") +
  axis(2, seq(0,100, by=20), col.axis="black", col = "black", las = 2, cex.axis=1.2) +
  axis(1, seq(-10,8, by=2), col.axis="black", col = "black", las = 1, cex.axis=1.2)




####### UPSETR plot --> Figure 3c ####### 

# Read file of up-regulated genes
up_genes = read.csv("upGenes_progressive_upsetR_TC_filter.csv", h=T,  na = "")

# Read it as data frame
df_up = as.data.frame(up_genes)

# View the list of named vectors to UpSetR converter
fromList(df_up)

# Plot upsetR
upset(fromList(df_up), line.size = 1, text.scale = 1.7, point.size=2.5, keep.order = T,
      sets = c("Up.T4.S.16.vs.T5.S.8", "Up.T3.S.25.vs.T4.S.16", "Up.T2.S.35.vs.T3.S.25" ,"Up.T1.S.45.vs.T2.S.35"),
      sets.bar.color=c("dodgerblue2" ,"goldenrod1", "darkorange2","firebrick"),
      main.bar.color = "black",
      queries = list(
        list(query = intersects, params = list("Up.T1.S.45.vs.T2.S.35"), color = "firebrick", active = T),
        list(query = intersects, params = list("Up.T3.S.25.vs.T4.S.16"), color = "goldenrod1", active = T),
        list(query = intersects, params = list("Up.T2.S.35.vs.T3.S.25"), color = "darkorange2", active = T),
        list(query = intersects, params = list("Up.T4.S.16.vs.T5.S.8"), color = "dodgerblue2", active = T)))

# Read file of down-regulated genes
down_genes = read.csv("downGenes_progressive_upsetR_TC_filter.csv", h=T)[,5:8]
df_down = as.data.frame(down_genes)
fromList(df_down)

upset(fromList(df_down), line.size = 1, text.scale = 1.7, point.size=2.5, keep.order = T,
      sets = c("Down.T4.S.16.vs.T5.S.8", "Down.T3.S.25.vs.T4.S.16", "Down.T2.S.35.vs.T3.S.25", "Down.T1.S.45.vs.T2.S.35"),
      sets.bar.color=c("turquoise3","dodgerblue2", "goldenrod1",  "darkorange2"),
      main.bar.color = "black", 
      queries = list(list(query = intersects, 
                          params = list("Down.T1.S.45.vs.T2.S.35"), color = "darkorange2", active = T),
                     list(query = intersects, params = list("Down.T2.S.35.vs.T3.S.25"), color = "goldenrod1", active = T),
                     list(query = intersects, params = list("Down.T3.S.25.vs.T4.S.16"), color = "dodgerblue2", active = T),
                     list(query = intersects, params = list("Down.T4.S.16.vs.T5.S.8"), color = "turquoise3", active = T)))


####### Scatter plot --> Figure 3d ####### 

# Read file of annotation
kegg_path_class = read.table("up_d_45_8_p_adj_KEGG.txt", h=T, sep="\t")
kog_path_class = read.table("up_d_45_8_p_adj_KOG.txt", h=T, sep="\t")

# Order labels of annotation
kegg_path_class = kegg_path_class[order(kegg_path_class$KEGG_pathways_class),]
kog_path_class = kog_path_class[order(kog_path_class$KOG_pathways_class),]

# Transform counts as integer
kegg_path_class$Gene_count = as.integer(kegg_path_class$Gene_count)
kog_path_class$Gene_count = as.integer(kog_path_class$Gene_count)

# Plot scatter
ggplot(kegg_path_class, aes(x=UP_DOWN, y=KEGG_pathways_class, color=UP_DOWN)) + 
  geom_point(aes(size=Gene_count)) +
  scale_y_discrete(limits = rev(levels(kegg_path_class$KEGG_pathways_class))) +
  scale_color_manual(values = c("firebrick3", "turquoise3"))

ggplot(kog_path_class, aes(x=UP_DOWN, y=KOG_pathways_class, color=UP_DOWN)) + 
  geom_point(aes(size=Gene_count)) +
  scale_y_discrete(limits = rev(levels(kog_path_class$KOG_pathways_class))) +
  scale_color_manual(values = c("firebrick3", "turquoise3"))



####### Canonical Analysis (CA) --> Figure 6c ####### 

# Read file of ice binding protein count in all samples
IBP_read = read.table("IBP_count.txt", row.names = 1, header=T)

# FILTER -- delete sample with 0 values 
Ophio_read = subset(Ophio_read, rowSums(Ophio_read)!=0)

# transpose
Ophio_read_t <- t(Ophio_read)

# metadata 
env = read.csv("Reads_per_gene/variables_reads.csv", row.names = 1, check.names = F)

# Verify if labels are ordered
env$Description == rownames(Ophio_read_t)

# If not ordered it
env<-env[rownames(Ophio_read_t),]

# Verify again
env$Description == rownames(Ophio_read_t)

# Standardize data
Ophio_read.hel = decostand(Ophio_read_t, method = "hellinger")

# Calcultate CA
ophio_ca = cca(Ophio_read.hel, dist="bray")

# ANOVA test
anova(ophio_ca, perm=999)

# Plot CA
plot(ophio_ca, type="p")
