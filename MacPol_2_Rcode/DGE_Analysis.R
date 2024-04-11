## DIFFERENTIAL EXPRESSION ANALYSIS

# Installing packages
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#BiocManager::install("mixOmics")
#BiocManager::install("HTSFilter")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("GO.db")
#BiocManager::install("biomaRt")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

library(limma)
library (edgeR)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(dplyr)
library(gplots)
library(heatmaply)
library(ggplot2)
library(ggrepel)

## Downloading required files (count- and metadata)-----------------------------

# Load in data
setwd("/Users/karinaovedal/Documents/DGE_Analysis/")

# Read in the raw count data
df <- read.delim("Count_DataFrame.csv", header=TRUE, sep=",")
#View(data)
head (data)

# Read the sample information
sampleInfo <- read.csv("Metadata.csv", header=TRUE, sep=",")
#View(sampleInfo)

## BCV analysis---------------------------------------------

# Create the DGEList
DGE <- DGEList(counts=df[,2:25], genes=df[,1],group=sampleInfo$Phenotype)

head (DGE) 
group <- factor(c("Phenotype"))

DGE$sample$group

# Filter to remove genes with counts less than 100
dim(DGE)
DGE.full <- DGE # keep the old one 
head(DGE$counts)
head(cpm(DGE))
apply(DGE$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(DGE)>100) >= 2 # Keeping only counts that are >100
DGE <- DGE[keep,]
dim(DGE) 

# Reset the library size after filtering
DGE$samples$lib.size <- colSums(DGE$counts)
DGE$samples

# Calculate normalization factors
DGE <- calcNormFactors(DGE)
DGE$samples

# Check the BCV (Biological coefficient of Variation) and dispersion of the samples

plotMDS(DGE, method="bcv", labels=sampleInfo$ID,col=as.numeric(DGE$samples$group))

y1 <- estimateCommonDisp(DGE, verbose=T)
y1 <- estimateTagwiseDisp(y1)
plotBCV(y1)

## GLM fitting-----------------------------

# Normalize genes for Gene Set Enrichment Analysis (GSEA) and save for analysis
y3<- as.matrix(y1)
write.csv(y3,"GSEA.csv", row.names = TRUE)

# Fit a generalized linear model (GLM) to the data to identify DEGs
design.mat <- model.matrix(~ 0 + DGE$samples$group)
colnames(design.mat) <- levels(DGE$samples$group)
design.mat


y2 <- estimateGLMCommonDisp(DGE,design.mat)
y2 <- estimateGLMTrendedDisp(y2,design.mat)
y2 <- estimateGLMTagwiseDisp(y2,design.mat)
names (y2)

# Create a scatter plot of the biological coefficient of variation (BCV) 
# against the average abundance of each gene
plotBCV(y2, col.common = "#e35a42", col.trend = "#5fb35b", col.tagwise = "#85513e")


## Prepare files for DEG analysis and heatmap----------------------------------- 

# Set up the contrast matrix
# Explicitly tell results which comparison to make

contrast_matrix <- makeContrasts(M1 = M11 - M01,
                                  M1P = M1P - M01,
                                  M2 = M22 - M02,
                                  M2P = M2P - M02,
                                  M2 = M23 - M13,
                                  M1M2P = M1M2P - M13,
                                  levels=design.mat)

contrast_matrix <- makeContrasts(M2 = M23 - M13,
                                 M1M2P = M1M2P - M13,
                                 levels=design.mat)
contrast_matrix
# ?makeContrasts

# make the file for creating the heatmap
design.mat

fit <- glmQLFit(y2, design.mat)

QLFT_all  <- glmQLFTest(fit, contrast=contrast_matrix)
topTags(QLFT_all)
QLFT_all<-topTags(QLFT_all, n=500, p.value =0.05) # Top five hundred genes
write.csv(QLFT_all,"heatmapfile.csv", row.names = FALSE)


## Creating the heatmap---------------------------------------------------------

x <- read.csv("heatmapfile.csv", header=TRUE, sep=",")

# prepare the data for the heatmap function
x = subset(x, select = -c(logCPM,F,PValue,FDR)) 
x.with.rownames <- data.frame(x[,-1], row.names=x[,1])
yh <- as.matrix (x.with.rownames)

# Row clustering (adjust here distance/linkage methods to what you need!)
hr <- hclust(as.dist(1-cor(t(yh), method="pearson")), method="average")

# Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(yh, method="pearson")), method="average")
dev.off()

# heat meap using pearson average, no scaling
# This part: labRow=rep("", nrow(yh)) removes the gene names from the heatmap
#heatmap.2(yh, 
#          Rowv=as.dendrogram(hr), 
#          Colv=as.dendrogram(hc), 
#          col=bluered, 
#          trace="none", 
#          margins=c(5,15), 
#          labRow=rep("", nrow(yh)), # Display row names
#          cexRow=0.9, # Increase this to make row labels (e.g., LogFC.M1) larger
#          cexCol=0.9 # Adjust if you also want to increase column label size
#)
# For switching, no need for the dendograms as there are only two columns.
heatmap.2(yh, 
          Rowv=as.dendrogram(hr), 
          Colv=FALSE, 
          col=bluered, 
          trace="none", 
          margins=c(5,15), 
          labRow=rep("", nrow(yh)),
          cexRow=0.7, # Increase this to make row labels (e.g., LogFC.M1) larger
          cexCol=0.7) # Adjust if you also want to increase column label size
         
## Other visualizations (static and interactive heatmaps)-------------------------

## Return matrix with row/column sorting as in heatmap- gives an error (subscript out of bounds)
#y[rev(hr$labels[hr$order]), hc$labels[hc$order]]

# Same heatmap from heatmap.2 as ggmap or heatmaply will produce 
heatmap.2(yh,trace="none", col=viridis(100), key=FALSE)

gradient_col <- ggplot2::scale_fill_gradient2(
  low = "blue", high = "red", 
  midpoint = 0, limits = c(-6, 6))

# Static version
ggheatmap(yh,scale_fill_gradient_fun = gradient_col, seriate = "mean", 
          row_dend_left = TRUE, k_col = 2, 
          k_row = 12)
# Interactive version
heatmaply(yh,scale_fill_gradient_fun = gradient_col, seriate = "mean", 
          row_dend_left = TRUE, k_col = 2, 
          k_row = 12)

# Save html file
dir.create("folder")
heatmaply(yh,scale_fill_gradient_fun = gradient_col, seriate = "mean", 
          row_dend_left = TRUE, k_col = 2, 
          k_row = 12, file = "folder/heatmaply_plot.html")
browseURL("folder/heatmaply_plot.html")


## DEG analysis using glmTreat--------------------------------------------------

# This section utilizes the glmTreat function to identify DEGs with 
# a specific log fold change (LFC) threshold, allowing for a more 
# stringent analysis of gene expression differences.
?glmTreat

# Using treat method for DGE- using a cut off of LFC 0.5, p= 0.05
design.mat
fit <- glmQLFit(y2, design.mat)

# "lfc= " specifies the minimum LFC threshold that
# that a gene must exceed to be considered significantly 
# differentially expressed. A higher LFC only keeps genes that are not only
# statistically significant but also biologically meaningful.
# Applying a stricter LFC criterion can help reduce the number of genes 
# identified as DE due to minor fluctuations that are statistically 
# significant but might not be biologically relevant = reducing false positives.

### Creating objects for the sample groups ###

#untreated_M1 <- glmTreat(fit, contrast=c(-1,0,1,0,0,0),lfc=log2(1.3)) 
treated_M1 <- glmTreat(fit, contrast=c(0, 0, -1, 0, 0, 1, 0, 0, 0),lfc=log2(1.3))
#untreated_M2 <- glmTreat(fit, contrast=c(0,0,-1,0,0,0),lfc=log2(1.3))
treated_M2 <- glmTreat(fit, contrast=c(0, 0, 0, 0, 0, 0, -1, 0, 1),lfc=log2(1.3))
treated_M2switch <- glmTreat(fit, contrast=c(0, 0, 0, 0, 1, 0, 0, -1, 0),lfc=log2(1.3))

### Creating objects for the sample groups - removing the lfc criteria ###

M1M2_noLFC <- glmTreat(fit, contrast=c(0, 0, 0, -1, 0, 0, 0, 1, 0))
M1M2P_noLFC <- glmTreat(fit, contrast=c(0, 0, 0, 0, 1, 0, 0, -1, 0))

# List all genes with LFC for M2_switch  
topTags_M1M2P <- topTags(treated_M2switch, n= 3571)
topTags(treated_M2switch)
topTags_M1M2P

# List all genes with LFC for treated_M1  
topTags_M1 <- topTags(treated_M1, n= 3571)
topTags(treated_M1)
topTags_M1


# List DEGs - treat analysis genes that have a log2 LFC absolute value 
# greater than 1.2 with a p-value cutoff of 0.05
# Benjamini-Hochberg method for adjusting p-values to control the 
# FDR, a common approach in multiple hypothesis testing to reduce 
# the chance of false positives.
DEG_M1 <- decideTestsDGE(treated_M1, adjust.method="BH", p.value = 0.05)
DEG_M2 <- decideTestsDGE(treated_M2, adjust.method="BH", p.value = 0.05)
DEG_M2switch <- decideTestsDGE(treated_M2switch, adjust.method="BH", p.value = 0.05)

DEG_M1M2 <- decideTestsDGE(M1M2_noLFC, adjust.method="BH", p.value = 0.05)
DEG_M1M2P <- decideTestsDGE(M1M2P_noLFC,adjust.method="BH", p.value = 0.05)

# Show number of DEGs
summary(DEG_M1)
summary(DEG_M2)
summary(DEG_M2switch)
summary(DEG_M1M2)
summary(DEG_M1M2P)

# List DEGs
DEGlist_M1 <- topTags(treated_M1, n= 3571, p.value =0.05)
DEGlist_M2 <- topTags(treated_M2, n= 3571, p.value =0.05)
DEGlist_M2switch <- topTags(treated_M2switch, n= 3571, p.value =0.05)

DEGlist_M1M2 <- topTags(M1M2_noLFC, n= 3571, p.value =0.05)
DEGlist_M1M2P <- topTags(M1M2P_noLFC, n= 3571, p.value =0.05)

# Write DEG lists to CSV files
write.csv(DEGlist_M1[["table"]], "/Users/karinaovedal/Documents/DGE_Analysis/DEG_M1.csv", row.names = FALSE)
write.csv(DEGlist_M2[["table"]], "/Users/karinaovedal/Documents/DGE_Analysis/DEG_M2.csv", row.names = FALSE)
write.csv(DEGlist_M2switch[["table"]], "/Users/karinaovedal/Documents/DGE_Analysis/DEG_M2switch.csv", row.names = FALSE)
write.csv(DEGlist_M1M2[["table"]], "/Users/karinaovedal/Documents/DGE_Analysis/DEG_M1M2.csv", row.names = FALSE)
write.csv(DEGlist_M1M2P[["table"]], "/Users/karinaovedal/Documents/DGE_Analysis/DEG_M1M2P.csv", row.names = FALSE)

## Venn Diagram of DEGs--------------------------------------------------------
dev.off()
# Combine the decision matrices for a Venn diagram analysis
DEG_Venn <- cbind(DEG_M1, DEG_M2, DEG_M2switch)
head(DEG_Venn)
write.csv(DEG_Venn, "/Users/karinaovedal/Documents/DGE_Analysis/DEG_Venn.csv", row.names = TRUE)

# Generate counts for Venn diagram
DEG_Venn_all <- vennCounts(DEG_Venn, include="both")

# Create the Venn diagram
vennDiagram(DEG_Venn_all, include="both", names=c("M1_DEG", "M2_DEG", "M2switch_DEG"), mar=rep(1,4), cex=c(1.5,1,0.7), lwd=1, circle.col=2, counts.col=NULL, show.include=NULL)

### Venn diagram of switch experiment: M2 vs M2+PGE2 ###
dev.off()
DEG_VennSwitch <- cbind(DEG_M1M2, DEG_M1M2P)

# Generate counts for Venn diagram
DEG_Venn_switch_all <- vennCounts(DEG_VennSwitch, include="both")

# Create the Venn diagram
vennDiagram(DEG_Venn_switch_all, include="both", names=c("M1M2_DEG", "M1M2P_DEG"), mar=rep(1,4), cex=c(1.5,0.7), lwd=1, circle.col=2, counts.col=NULL, show.include=NULL)


# Extract gene names from DEG lists (second column)
rownames(DEGlist_M1[["table"]]) <- DEGlist_M1[["table"]][,1]  
rownames(DEGlist_M2[["table"]]) <- DEGlist_M2[["table"]][,1]
rownames(DEGlist_M2switch[["table"]]) <- DEGlist_M2switch[["table"]][,1]
rownames(DEGlist_M1M2[["table"]]) <- DEGlist_M1M2[["table"]][,1]

genes_M1 <- rownames(DEGlist_M1[["table"]])
genes_M2 <- rownames(DEGlist_M2[["table"]])
genes_M2switch <- rownames(DEGlist_M2switch[["table"]])
genes_M2untr <- rownames(DEGlist_M1M2[["table"]])

# Find common genes between M1 and M2 macrophages
common_M1_M2 <- intersect(genes_M1, genes_M2)

# Find common genes between M1 and M2 switch macrophages
common_M1_M2switch <- intersect(genes_M1, genes_M2switch)

# Find common genes between M2 and M2 switch macrophages
common_M2_M2switch <- intersect(genes_M2, genes_M2switch)

# Find common genes between untreated versus treated M2 in switch experiment
common_M2_M2 <- intersect(genes_M2untr, genes_M2switch)

# Find genes common to all groups
common_all <- Reduce(intersect, list(genes_M1, genes_M2, genes_M2switch))

# Print the results
cat("Common genes between M1 and M2 macrophages:\n")
print(common_M1_M2)
cat("\nCommon genes between M1 and M2 switch macrophages:\n")
print(common_M1_M2switch)
cat("\nCommon genes between untreated versus treated M2 in switch experiment:\n")
print(common_M2_M2)
cat("\nGenes common to all groups:\n")
print(common_all)

## VOLCANO PLOT (M2 switch)---------------------------------------------------------------
# If you want the first heatmap to run, the make sure to not run the code here.
dev.off()

deg_results <- topTags(treated_M2switch, n=Inf)  # n=Inf to get all results
volcano_data <- with(deg_results$table, data.frame(logFC=logFC, pValue=PValue, Gene=genes)) 

# Marking the top 41 significant genes based on p-value
volcano_data$topGene <- FALSE
topGenes <- volcano_data[volcano_data$pValue < 0.05,]
topGenes <- topGenes[order(topGenes$pValue),][1:20,]
volcano_data$topGene[rownames(volcano_data) %in% rownames(topGenes)] <- TRUE
head(volcano_data)

# Create the volcano plot
ggplot(volcano_data, aes(x=logFC, y=-log10(pValue), label=Gene)) +
  geom_point(aes(color=topGene), alpha=0.5) + 
  scale_color_manual(values=c(`FALSE`="#85513e", `TRUE`="#29ba14")) + 
  geom_text_repel(data=subset(volcano_data, topGene==TRUE), 
                  aes(label=Gene), size=4.5,
                  box.padding = 0.35, point.padding = 0.6,
                  max.overlaps = Inf) +
  labs(x="Log2 Fold Change", y="-Log10 P-value", title="RE-POLARIZATION: \nUntreated vs. PGE2 Treated M2 macrophages") +
  theme_minimal() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, face="bold", size=20, color="#85513e"),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#e35a42") +
  scale_y_continuous(breaks=seq(0, max(-log10(volcano_data$pValue)), by=1))

## VOLCANO PLOT (M1)---------------------------------------------------------------
# If you want the first heatmap to run, the make sure to not run the code here.
dev.off()

deg_results <- topTags(treated_M1, n=Inf)  # n=Inf to get all results
volcano_data <- with(deg_results$table, data.frame(logFC=logFC, pValue=PValue, Gene=genes)) 

# Marking the top 20 significant genes based on p-value
volcano_data$topGene <- FALSE
topGenes <- volcano_data[volcano_data$pValue < 0.05,]
topGenes <- topGenes[order(topGenes$pValue),][1:20,]
volcano_data$topGene[rownames(volcano_data) %in% rownames(topGenes)] <- TRUE

# Create the volcano plot for M1 + PGE2
ggplot(volcano_data, aes(x=logFC, y=-log10(pValue), label=Gene)) +
  geom_point(aes(color=topGene), alpha=0.5) + 
  scale_color_manual(values=c("FALSE"="#85513e", "TRUE"="#29ba14")) +
  geom_text_repel(data=subset(volcano_data, topGene==TRUE), 
                  aes(label=Gene), size=4.5,
                  box.padding = 0.35, point.padding = 0.5,
                  max.overlaps = Inf) +
  labs(x="Log2 Fold Change", y="-Log10 P-value", title="Untreated vs. PGE2 Treated M1") + 
  theme_minimal() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, face="bold", size=20, color="#85513e"),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#e35a42") +
  scale_y_continuous(breaks=seq(0, max(-log10(volcano_data$pValue)), by=1)) +
  xlim(c(-6,6))

## VOLCANO PLOT (M2)---------------------------------------------------------------
# If you want the first heatmap to run, the make sure to not run the downstream code here.
dev.off()

deg_results <- topTags(treated_M2, n=Inf)  # n=Inf to get all results
volcano_data <- with(deg_results$table, data.frame(logFC=logFC, pValue=PValue, Gene=genes)) 

# Marking the top 20 significant genes based on p-value
volcano_data$topGene <- FALSE
topGenes <- volcano_data[volcano_data$pValue < 0.05,]
topGenes <- topGenes[order(topGenes$pValue),][1:20,]
volcano_data$topGene[rownames(volcano_data) %in% rownames(topGenes)] <- TRUE

# Create the volcano plot for M2 + PGE2
ggplot(volcano_data, aes(x=logFC, y=-log10(pValue), label=Gene)) +
  geom_point(aes(color=topGene), alpha=0.5) + 
  scale_color_manual(values=c(`FALSE`="#85513e", `TRUE`="#29ba14")) +
  geom_text_repel(data=subset(volcano_data, topGene==TRUE), 
                  aes(label=Gene), size=4.5,
                  box.padding = 0.35, point.padding = 0.5,
                  max.overlaps = Inf) +
  labs(x="Log2 Fold Change", y="-Log10 P-value", title="Untreated vs. PGE2 Treated M2") + 
  theme_minimal() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, face="bold", size= 20, color="#85513e"),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#e35a42") +
  scale_y_continuous(breaks=seq(0, max(-log10(volcano_data$pValue)), by=1)) +
  xlim(c(-10,10))
    


