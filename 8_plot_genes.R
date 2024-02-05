### STEP 6: Plot significant DE genes using multiple visualizations ###

load("~/opt/r/rnaseq-riboseq/RNA-RIBO_DESeqResults.RData")

# Load Libraries
library(tidyverse)
library(AnnotationDbi)
library(ggrepel)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(cowplot)
library(EnhancedVolcano)
library(org.Hs.eg.db)

# Single-gene Expression --------------------------------------------------
source("singlegeneExpr.R")
symbol <- "VEGFA"
gene_id <- mapIds(org.Hs.eg.db,
                  keys = symbol,
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  multiVals = "first")
exprs <- singlegeneExpr(gene_id, rna_dds, rpf_dds, symbol)

plot_grid(
  exprs[["rna"]],
  exprs[["rpf"]],
  align = 'h',
  ncol = 2,
  labels = "AUTO",
  vjust = 1
)

# Plot counts of top 20 DE Genes ---------------------------------------------------------
# Order results by padj values
top20genes_rna <- resdf_rna %>% 
  arrange(padj) %>%   #Arrange rows by padj values
  pull(entrez) %>% 	#Extract character vector of ordered genes
  head(n=20)		      #Extract the first 20 genes

top20genes_rpf <- resdf_rpf %>% 
  arrange(padj) %>%
  pull(entrez) %>%
  head(n=20)

# Fetch normalized counts for top 20 significant genes
top20norm_rna <- normCounts_rna %>%
  dplyr::filter(gene_id %in% top20genes_rna)
# append gene symbol column to end of dataframe
top20norm_rna <- left_join(x = top20norm_rna, y = dplyr::select(resdf_rna, entrez, symbol), by = c("gene_id" = "entrez"))
# adjust column names in normalized matrix to sure proper joining after data manipulation below
colnames(top20norm_rna)[2:9] <- pull(sampleTable.rna, sample) 

top20norm_rpf <- normCounts_rpf %>%
  dplyr::filter(gene_id %in% top20genes_rpf)
top20norm_rpf <- left_join(x = top20norm_rpf, y = dplyr::select(resdf_rpf, entrez, symbol), by = c("gene_id" = "entrez"))
colnames(top20norm_rpf)[2:9] <- pull(sampleTable.rpf, sample)

# Gather the columns to have normalized counts to a single column
top20norm_rna <- top20norm_rna %>%
  gather(colnames(top20norm_rna)[2:9], key = "samplename", value = "normalized_counts")

top20norm_rpf <- top20norm_rpf %>%
  gather(colnames(top20norm_rpf)[2:9], key = "samplename", value = "normalized_counts")

# Check the column header and sample names in the "gathered" data frame
View(top20norm_rna)
View(sampleTable.rna)

top20norm_rna <- left_join(top20norm_rna, sampleTable.rna, by = c("samplename" = "sample"))
top20norm_rpf <- left_join(top20norm_rpf, sampleTable.rpf, by = c("samplename" = "sample"))

# Plot top 20 DE genes using ggplot2 
ggplot(top20norm_rna) +
  geom_point(aes(x = symbol, y = normalized_counts, color = treatment)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  labs(color = "Treatment") +
  ggtitle("Top 20 DE Genes in RNA-seq") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(top20norm_rpf) +
  geom_point(aes(x = symbol, y = normalized_counts, color = treatment)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  labs(color = "Condition") +
  ggtitle("Top 20 DE Genes in RIBO-seq") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

# Volcano Plots ------------------------------------------------------------
volc_rna <- EnhancedVolcano(resKeepOrdered_rna, # Width: 800, Height: 600
                            lab = resKeepOrdered_rna$symbol,
                            x = 'log2FoldChange',
                            y = 'padj',
                            #title = "Volcano plot of RNA-seq",
                            #subtitle = "Hypoxia versus Normoxia",
                            pCutoff = 1e-05,
                            FCcutoff = 1,
                            legendPosition = "right")

volc_rpf <- EnhancedVolcano(resKeepOrdered_rpf,
                            lab = resKeepOrdered_rpf$symbol,
                            x = 'log2FoldChange',
                            y = 'padj',
                            #title = "Volcano plot of RIBO-seq",
                            #subtitle = "Hypoxia versus Normoxia",
                            pCutoff = 1e-05,
                            FCcutoff = 1,
                            legendPosition = "right")

# Extract the legend from one of the plots
legend <- get_legend(
  volc_rna + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# Arrange the plots in a single row
prow <- plot_grid(
  volc_rna + theme(legend.position = "none"),
  volc_rpf + theme(legend.position = "none"),
  align = 'vh',
  labels = "AUTO",
  nrow = 1
)

# Add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(prow, legend, rel_widths = c(3, .5))

# Heatmaps -----------------------------------------------------------------
# Extract normalized expression for all significant genes (padj)
padj.cutoff <- 0.10 # set threshold for FDR
sigdf_rna <- resdf_rna %>%
  filter(padj < padj.cutoff)
sigdf_rpf <- resdf_rpf %>%
  filter(padj < padj.cutoff)
signorm_rna <- normCounts_rna %>% 
  dplyr::filter(gene_id %in% sigdf_rna$entrez)
signorm_rpf <- normCounts_rpf %>% 
  dplyr::filter(gene_id %in% sigdf_rpf$entrez)

# Isolate all significant genes with LFC one sigma above the mean
sd_rna <- sd(resKeepOrdered_rna$log2FoldChange, na.rm = TRUE)
mean_rna <- mean(resKeepOrdered_rna$log2FoldChange, na.rm = TRUE)
poscutoff_rna <- mean_rna + (1.96 * sd_rna)
num_rna <- resdf_rna %>% filter(log2FoldChange > poscutoff_rna, padj < padj.cutoff)

# Isolate top 50 significant genes with most positive and most negative LFC,
# after sorting for the top 500 DE genes
pos50_sig.rna <- resdf_rna %>%
  filter(padj < padj.cutoff) %>%
  arrange(padj) %>%
  head(n=500) %>%
  arrange(desc(log2FoldChange)) %>%
  drop_na(symbol) %>%
  head(n=50)
neg50_sig.rna <- resdf_rna %>%
  filter(padj < padj.cutoff) %>%
  arrange(padj) %>%
  head(n=500) %>%
  arrange(log2FoldChange) %>%
  drop_na(symbol) %>%
  head(n=50)
all50_sig.rna <- resdf_rna %>%
  filter(padj < padj.cutoff) %>%
  arrange(padj) %>%
  head(n=500) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  drop_na(symbol) %>%
  head(n=50)

pos50_sig.rpf <- resdf_rpf %>%
  filter(padj < padj.cutoff) %>%
  arrange(padj) %>%
  head(n=500) %>%
  arrange(desc(log2FoldChange)) %>%
  drop_na(symbol) %>%
  head(n=50)
neg50_sig.rpf <- resdf_rpf %>%
  filter(padj < padj.cutoff) %>%
  arrange(padj) %>%
  head(n=500) %>%
  arrange(desc(log2FoldChange)) %>%
  drop_na(symbol) %>%
  tail(n=50)
all50_sig.rpf <- resdf_rpf %>%
  filter(padj < padj.cutoff) %>%
  arrange(padj) %>%
  head(n=500) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  drop_na(symbol) %>%
  head(n=50)

# Extract normalized expression for top/bottom 50 LFC genes from the samples
pos50_norm.rna <- normCounts_rna %>% 
  dplyr::filter(gene_id %in% pos50_sig.rna$entrez)
neg50_norm.rna <- normCounts_rna %>% 
  dplyr::filter(gene_id %in% neg50_sig.rna$entrez)
all50_norm.rna <- normCounts_rna %>% 
  dplyr::filter(gene_id %in% all50_sig.rna$entrez)


pos50_norm.rpf <- normCounts_rpf %>% 
  dplyr::filter(gene_id %in% pos50_sig.rpf$entrez)
neg50_norm.rpf <- normCounts_rpf %>% 
  dplyr::filter(gene_id %in% neg50_sig.rpf$entrez)
all50_norm.rpf <- normCounts_rpf %>% 
  dplyr::filter(gene_id %in% all50_sig.rpf$entrez)

# Set heatmap parameters
heat_colors <- brewer.pal(6, "YlOrRd")
heat_annot.rna <- sampleTable.rna %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select("Cell Line" = cell_line, "Treatment" = treatment)
rownames(heat_annot.rna) <- colnames(all50_norm.rna)[2:9]
heat_annot.rpf <- sampleTable.rpf %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select("Cell Line" = cell_line, "Treatment" = treatment)
rownames(heat_annot.rpf) <- colnames(all50_norm.rpf)[2:9]

all50_norm.rna$symbol <- mapIds(org.Hs.eg.db,
                         keys = all50_norm.rna$gene_id,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
all50_norm.rpf$symbol <- mapIds(org.Hs.eg.db,
                         keys = all50_norm.rpf$gene_id,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")

heat_rownames_all.rna <- all50_norm.rna$symbol
heat_rownames_all.rpf <- all50_norm.rpf$symbol

samp_annot.rna <- heat_annot.rna %>% arrange(Treatment)
samp_ordered.rna <- all50_norm.rna[rownames(samp_annot.rna)]

heat.rna <- pheatmap(samp_ordered.rna, # Width: 600, Height: 1000
                     color = heat_colors,
                     cluster_rows = T,
                     labels_row = heat_rownames_all.rna,
                     show_rownames = T,
                     show_colnames = F,
                     annotation = samp_annot.rna,
                     cluster_cols = F,
                     border_color = NA,
                     fontsize = 10,
                     scale = "row",
                     fontsize_row = 6,
                     height = 20)

samp_annot.rpf <- heat_annot.rpf %>% arrange(Treatment)
samp_ordered.rpf <- all50_norm.rpf[rownames(samp_annot.rpf)]

heat.rpf <- pheatmap(samp_ordered.rpf,
                     color = heat_colors,
                     cluster_rows = T,
                     labels_row = heat_rownames_all.rpf,
                     show_rownames = T,
                     show_colnames = F,
                     annotation = samp_annot.rpf,
                     cluster_cols = F,
                     border_color = NA,
                     fontsize = 10,
                     scale = "row",
                     fontsize_row = 6,
                     height = 20)

# Arrange the plots in a single row
prow <- plot_grid(
  heat.rna[[4]],
  heat.rpf[[4]],
  align = 'vh',
  labels = "AUTO",
  nrow = 1
)

# Add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(prow)