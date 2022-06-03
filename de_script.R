### Gene-level differential expression analysis using DESeq2

# load libraries
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(ggplot2)
library(apeglm)
library(ggrepel)

# load data
data <- read.table("data/RNA_featurecounts.Rmatrix.txt", header = T, row.names = 1)
meta <- read.table("meta/count_tab_sample.txt", header = T, row.names = 1) %>% 
  slice(1:8) %>%
  mutate_at(vars(cellType, condition), as.factor)

# plot count distribution
ggplot(data) +
  geom_histogram(aes(x = S1), stat = "bin", bins = 200) +
  xlim(-5, 500) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# modeling count data of GSC215 hypoxic replicates
mean_counts <- apply(data[, 3:4], 1, mean)
variance_counts <- apply(data[, 3:4], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) +
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()

# match metadata and counts data
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta,
                              design = ~ cellType + condition)
## Generate normalized counts

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/RNA_normalized_counts.txt",
            sep="\t",
            quote=F,
            col.names=NA)

## Quality Assessment using DESeq2

# transform counts for data visualization
rld <- rlog(dds, blind = TRUE)

# plot pca
plotPCA(rld, intgroup="condition")
plotPCA(rld, intgroup="cellType")

# extract the rlog matrix from the object
rld_mat <- assay(rld)

# compute pairwise correlation values
rld_cor <- cor(rld_mat)
head(rld_cor)
head(meta)

# plot heatmap
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation = meta,
         color = heat.colors,
         border_color=NA,
         fontsize = 10, 
         fontsize_row = 10,
         height=20)

# run analysis
dds <- DESeq(dds)

# plot dispersion estimates
plotDispEsts(dds)

## Exploring results

results(dds)
contrast_hn <- c("condition", "Hypoxic", "Normoxic")
res_tableHN <- results(dds, contrast = contrast_hn, alpha = 0.05)
contrast_ct <- c("cellType", "315", "216")
res_tableCT <- results(dds, contrast = contrast_ct, alpha = 0.05)

# Check what type of object is returned
class(res_tableHN)
class(res_tableCT)

# What is stored in results?
res_tableHN %>% 
  data.frame() %>% 
  View()

# Get information on each column in results
mcols(res_tableHN, use.names=T)
mcols(res_tableCT, use.names=T)

## Fold Change

# Save the unshrunken results to compare
res_tableHN_unshrunken <- res_tableHN
res_tableCT_unshrunken <- res_tableCT

# Apply fold change shrinkage
res_tableHN <- lfcShrink(dds,
                       coef="condition_Normoxic_vs_Hypoxic",
                       type="apeglm")
res_tableCT <- lfcShrink(dds,
                           coef="cellType_315_vs_216",
                           type="apeglm")

# MA plot using unshrunken fold changes
plotMA(res_table_unshrunken, ylim=c(-2,2))

# MA plot using shrunken fold changes
plotMA(res_tableHN, ylim=c(-2,2))
plotMA(res_tableCT, ylim=c(-2,2))

## Summarize results
summary(res_tableHN, alpha = 0.05)
summary(res_tableCT, alpha = 0.05)

# Set thresholds
padj.cutoff <- 0.05

# Create a tibble of results
res_tableHN_tb <- res_tableHN %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_tableCT_tb <- res_tableCT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the tibble to keep only significant genes
sigHN <- res_tableHN_tb %>%
  filter(padj < padj.cutoff)
sigCT <- res_tableCT_tb %>%
  filter(padj < padj.cutoff)

## Visualizing the results

RNA_meta <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

# First convert normalized_counts to a data frame and transfer the row names to a new column called "gene"
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

# Next, merge together (ensembl IDs) the normalized counts data frame with a subset of the annotations in the tx2gene data frame (only the columns for ensembl gene IDs and gene symbols)
grch38annot <- annotations_ahb %>% 
  dplyr::select(gene_id_version, gene_name) %>% 
  dplyr::distinct()

# This will bring in a column of gene symbols
normalized_counts <- merge(normalized_counts, grch38annot, by.x="gene", by.y="gene_id_version")

# Now create a tibble for the normalized counts
normalized_counts <- normalized_counts %>%
  as_tibble()

normalized_counts 

## Plotting significant DE genes

# Find the Ensembl ID of HIF
grch38annot[grch38annot$gene_name == "HIF1A", "gene_id_version"]

# Plot expression for single gene
plotCounts(dds, gene="ENSG00000100644.17", intgroup="condition") 

# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="ENSG00000100644.17", intgroup="condition", returnData=TRUE)

# What is the data output of plotCounts()?
d %>% View()

# Plot the HIF1A normalized counts, using the samplenames (rownames(d) as labels)
ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("HIF1A") +
  theme(plot.title = element_text(hjust = 0.5))

## Plot top 20 genes

# Order results by padj values
top20_sigHN_genes <- res_tableHN_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20)		#Extract the first 20 genes

# normalized counts for top 20 significant genes
top20_sigHN_norm <- normalized_counts %>%
  filter(gene %in% top20_sigHN_genes)

# Gathering the columns to have normalized counts to a single column
gathered_top20_sigHN <- top20_sigHN_norm %>%
  gather(colnames(top20_sigHN_norm)[2:9], key = "samplename", value = "normalized_counts")

# check the column header in the "gathered" data frame
View(gathered_top20_sigHN)

gathered_top20_sigHN <- inner_join(RNA_meta, gathered_top20_sigHN)

# plot using ggplot2
ggplot(gathered_top20_sigHN) +
  geom_point(aes(x = gene_name, y = normalized_counts, color = condition)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

## Heatmap

# Extract normalized expression for significant genes from the OE and control samples (2:4 and 7:9)
norm_HNsig <- normalized_counts %>% 
  filter(gene %in% sigHN$gene) 

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_HNsig[2:9], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

## Volcano Plot

# Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction

res_tableHN_tb <- res_tableHN_tb %>% 
  mutate(threshold_HN = padj < 0.05 & abs(log2FoldChange) >= 0.58)

# plot
ggplot(res_tableHN_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_HN)) +
  ggtitle("Hypoxia vs. Normoxia") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

# Add all the gene symbols as a column from the grch38 table using bind_cols()
res_tableHN_tb <- bind_cols(res_tableHN_tb, symbol=grch38annot$gene_name[match(res_tableHN_tb$gene, grch38annot$gene_id_version)])

# Create an empty column to indicate which genes to label
res_tableHN_tb <- res_tableHN_tb %>% mutate(genelabels = "")

## Sort by padj values 
res_tableHN_tb <- res_tableHN_tb %>% arrange(padj)

## Populate the genelabels column with contents of the gene symbols column for the first 10 rows, i.e. the top 10 most significantly expressed genes
res_tableHN_tb$genelabels[1:10] <- as.character(res_tableHN_tb$symbol[1:10])

View(res_tableHN_tb)

ggplot(res_tableHN_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_HN)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Hypoxia vs. Normixia DE") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


