## This script performs pairwise differential gene expression analysis between pre- and mid- and between pre- and post-infection samples. Results from this script are presented under " Expression dynamics of ADAR genes and isofroms" in results.

## The script generates supplementary tables : 1A, 1B, 1D, and 1E

#load library
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Set the working directory
directory <- "F:\\Manuscript_1\\DESeq2"
setwd(directory)

# input count data
data  <- as.matrix(read.csv("gene_count_matrix.csv", header=T, row.names=1, sep=","))
head(data)

# create experiment labels 
colData <- read.csv("colData.csv", header=T, row.names=1, sep=",")
head(colData)

# making sure that the row names in col data match the column names in count data 
all(colnames(data) %in% rownames(colData))

# making sure that the row names in col data are in the same order as the count data 
all(colnames(data) == rownames(colData))

# create DESeq input matrix  
dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData = colData, 
                              design = ~infection)
dds

# Run DESeq 
dds <- DESeq(dds)
res <- results(dds)
res

## results for pre vs mid 
outputPrefix <- "Aiswarya_pre vs mid _DESeq2__"

## pairwise comparisons - Pre vs Mid
results_Pre_vs_mid <- results(dds, contrast = c("infection","Mid","Controls" ))
results_Pre_vs_mid
write.csv(results_Pre_vs_mid, file = paste0(outputPrefix, "-gene_counts_pre_vs_mid_before_filtering.csv"))

# applying filters for padj and FC 
ressubset_Pre_vs_mid_0.5 = subset(results_Pre_vs_mid,  abs(results_Pre_vs_mid$log2FoldChange) >= 0.58 & padj<0.05) 
ressubset_Pre_vs_mid_0.5 <- ressubset_Pre_vs_mid_0.5[order(ressubset_Pre_vs_mid_0.5$padj),]
write.csv(ressubset_Pre_vs_mid_0.5, file = paste0(outputPrefix, "-Gene_counts_pre_vs_mid_after_filters.csv"))

## results for pre vs post
outputPrefix <- "Aiswarya_pre vs post _DESeq2__"

## pairwise comparisons - Pre vs post
results_Pre_vs_post <- results(dds, contrast = c("infection","Post","Controls"))
results_Pre_vs_post
write.csv(results_Pre_vs_post, file = paste0(outputPrefix, "-gene_counts_pre_vs_post_before_filtering.csv"))

# applying filters for padj and FC 
ressubset_Pre_vs_post_0.5 = subset(results_Pre_vs_post,  abs(results_Pre_vs_post$log2FoldChange) >= 0.58 & padj<0.05) 
ressubset_Pre_vs_post_0.5 <- ressubset_Pre_vs_post_0.5[order(ressubset_Pre_vs_post_0.5$padj),]
write.csv(ressubset_Pre_vs_post_0.5, file = paste0(outputPrefix, "-Gene_counts_pre_vs_post_after_filters.csv"))


#------------------------------------------------------------------------------------------------------------------

## Transcript count matrix

# input count data
data  <- as.matrix(read.csv("transcript_count_matrix.csv", header=T, row.names=1, sep=","))
head(data)

# create experiment labels 
colData <- read.csv("colData.csv", header=T, row.names=1, sep=",")
head(colData)

# making sure that the row names in col data match the column names in count data 
all(colnames(data) %in% rownames(colData))

# making sure that the row names in col data are in the same order as the count data 
all(colnames(data) == rownames(colData))

# create DESeq input matrix  
dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData = colData, 
                              design = ~infection)
dds

# Run DESeq 
dds <- DESeq(dds)
res <- results(dds)
res

## results for pre vs mid 
outputPrefix <- "Aiswarya_pre vs mid _DESeq2__"

## pairwise comparisons - Pre vs Mid
results_Pre_vs_mid <- results(dds, contrast = c("infection","Mid","Controls" ))
results_Pre_vs_mid
write.csv(results_Pre_vs_mid, file = paste0(outputPrefix, "-transcript_counts_pre_vs_mid_before_filtering.csv"))


## results for pre vs post
outputPrefix <- "Aiswarya_pre vs post _DESeq2__"

## pairwise comparisons - Pre vs post
results_Pre_vs_post <- results(dds, contrast = c("infection","Post","Controls"))
results_Pre_vs_post
write.csv(results_Pre_vs_post, file = paste0(outputPrefix, "-transcript_counts_pre_vs_post_before_filtering.csv"))


outputPrefix <- "Aiswarya_DESeq2_transcript_"
normalized_transcript_expression <- as.data.frame(counts(dds,normalized =TRUE)) # this will include samples from all the conditions
normalized_transcript_expression # This is a dataframe with all the samples and normlalized counts
write.csv(normalized_transcript_expression, file = paste0(outputPrefix, "-normalized_counts.csv")) # use to make ADARp150vs adarp110 graphs

# plotting
# Only ADARp150 (ENST00000368474) and ADARp110 (ADARp150368471) normalized values were extracted from the above table and saved as a .csv file 

## ADAR1 plots

ADAR_isoform <- read.csv("ADAR_isoform.csv")
head(ADAR_isoform)


# Reorder infection stage 
ADAR_isoform$Infection_stage <- factor(ADAR_isoform$Infection_stage, levels = c("Controls", "Mid", "Post"))

p1 <- ggplot(ADAR_isoform, aes(x = Infection_stage, y = ADARp110, color = Infection_stage)) +
  geom_boxplot(fill = "white", alpha = 0.5) +  
  geom_jitter(width = 0.1, alpha = 0.9) +
  labs(x = "Stages of SARS-CoV-2 infection", y = "ADARp110 expression normalized") +
  scale_color_manual(values = c("darkmagenta", "red", "blue")) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 8)) +
  scale_x_discrete(labels = c("Controls" = "Pre-infection", "Mid" = "Mid-infection", "Post" = "Post-infection"))

# adding stats (paired t-test)
ADARp110 <- p1 + theme(legend.position = "none")
ADARp110 <- ADARp110 + stat_compare_means(comparisons = list(c("Controls", "Post"), c("Controls", "Mid")), method = "t.test", paired = TRUE)


p2 <- ggplot(ADAR_isoform, aes(x = Infection_stage, y = ADARp150, color = Infection_stage)) +
  geom_boxplot(fill = "white", alpha = 0.5) +  
  geom_jitter(width = 0.1, alpha = 0.9) +
  labs(x = "Stages of SARS-CoV-2 infection", y = "ADARp150 expression normalized") +
  scale_color_manual(values = c("darkmagenta", "red", "blue")) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 8)) +
  scale_x_discrete(labels = c("Controls" = "Pre-infection", "Mid" = "Mid-infection", "Post" = "Post-infection"))

# adding stats (paired t-test)
ADARp150 <- p2 + theme(legend.position = "none")
ADARp150 <- ADARp150 + stat_compare_means(comparisons = list(c("Controls", "Post"), c("Controls", "Mid")), method = "t.test", paired = TRUE)


# Arrange and label in grid
ggarrange(ADARp110, ADARp150, ncol = 2, labels= c("A", "B"))


