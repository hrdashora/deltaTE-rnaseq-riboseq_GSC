### STEP 5: Plotting significant DE genes ###
load("~/opt/r/rnaseq-riboseq/RNA-RIBO_DESeqResults.RData")

# Load Libraries
library(tidyr)
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
symbol <- "ALDOA"
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

tcga <- readr::read_csv("data/TCGA_GBM_CDKN2C_Exprs")
  
ggplot(tcga, aes(x = Histology, y = mRNA)) +
  geom_boxplot(outlier.alpha = 0.1) +
  ggpubr::stat_compare_means(mapping = aes(x = Histology, y = mRNA),
                             hide.ns = T,
                             method = "wilcox.test",
                             paired = F) +
  labs(y = "mRNA expression (log2)") +
  theme_classic()

cgga <- readr::read_csv("data/CGGA_CDKN2C_Exprs")

ggplot(cgga, aes(x = Histology, y = mRNA)) +
  geom_boxplot(outlier.alpha = 0.1) +
  ggpubr::stat_compare_means(mapping = aes(x = Histology, y = mRNA),
                             hide.ns = T,
                             paired = F) +
  labs(y = "mRNA expression (log2)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

# Visualize the expression profile
ggpubr::ggboxplot(cgga, x = "Histology", y = "mRNA", color = "Histology", legend = "none") +
  ggpubr::rotate_x_text(angle = 45) +
  geom_hline(yintercept = mean(cgga$mRNA), linetype = 2)+ # Add horizontal line at base mean
  ggpubr::stat_compare_means(label.y = 10)+ # Add global p-value
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.")        

# Top 20 DE Genes ---------------------------------------------------------
# Order results by padj values
top20_genes.rna <- resdf.rna %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene_id) %>% 		#Extract character vector of ordered genes
  head(n=20)		#Extract the first 20 genes

top20_genes.rpf <- resdf.rpf %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene_id) %>% 		#Extract character vector of ordered genes
  head(n=20)		#Extract the first 20 genes

# Normalized counts for top 20 significant genes
top20_norm.rna <- normCounts.rna %>%
  dplyr::filter(gene_id %in% top20_genes.rna)

top20_norm.rpf <- normCounts.rpf %>%
  dplyr::filter(gene_id %in% top20_genes.rna)

# Gather the columns to have normalized counts to a single column
top20_norm.rna <- top20_norm.rna %>%
  gather(colnames(top20_norm.rna)[2:9], key = "samplename", value = "normalized_counts")

top20_norm.rpf <- top20_norm.rpf %>%
  gather(colnames(top20_norm.rpf)[2:9], key = "samplename", value = "normalized_counts")

# Check the column header in the "gathered" data frame
View(top20_norm.rna)

top20_norm.rna <- inner_join(top20_norm.rna, sampleTable.rna, by = c("samplename" = "sample"))
top20_norm.rpf <- inner_join(top20_norm.rpf, sampleTable.rpf, by = c("samplename" = "sample"))

# Plot top 20 DE genes using ggplot2 
ggplot(top20_norm.rna) +
  geom_point(aes(x = SYMBOL, y = normalized_counts, color = treatment)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  labs(color = "Treatment") +
  ggtitle("Top 20 Significant DE Genes in RNA-seq") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(top20_norm.rpf) +
  geom_point(aes(x = SYMBOL, y = normalized_counts, color = treatment)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  labs(color = "Condition") +
  ggtitle("Top 20 Significant DE Genes in RIBO-seq") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

# Volcano Plot ------------------------------------------------------------

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

# Heatmap -----------------------------------------------------------------
# # Unannotated heatmap of 30 significant genes
# vsdMat.rna <- vst(ddsMat.rna)
# mat.rna <- assay(vsdMat.rna)[head(order(resdf.rna$padj),30),]
# mat.rna <- mat.rna - rowMeans(mat.rna)
# df <- as.data.frame(colData(vsdMat.rna)[,c("cell_line","treatment")])
# pheatmap(mat.rna, annotation_col=df)
# 
# vsdMat.rpf <- vst(ddsMat.rpf)
# mat.rpf <- assay(vsdMat.rpf)[head(order(resdf.rpf$padj),30),]
# mat.rpf <- mat.rpf - rowMeans(mat.rpf)
# df <- as.data.frame(colData(vsdMat.rpf)[,c("cell_line","treatment")])
# pheatmap(mat.rpf, annotation_col=df)

# Extract normalized expression for all significant genes (padj)
# Set thresholds
padj.cutoff <- 0.05
sigdf_rna <- resdf_rna %>%
  filter(padj < padj.cutoff)
sigdf_rpf <- resdf_rpf %>%
  filter(padj < padj.cutoff)
signorm_rna <- normCounts_rna %>% 
  dplyr::filter(gene_id %in% sigdf_rna$entrez) 
signorm_rpf <- normCounts_rpf %>% 
  dplyr::filter(gene_id %in% sigdf_rpf$entrez)

# Isolate top 50 significant genes with most positive and most negative LFC,
# after sorting for the top 500 DE genes
sd_rna <- sd(resKeepOrdered_rna$log2FoldChange, na.rm = TRUE)
mean_rna <- mean(resKeepOrdered_rna$log2FoldChange, na.rm = TRUE)
poscutoff_rna <- mean_rna + (1.96 * sd_rna)
num_rna <- resdf_rna %>% filter(log2FoldChange > poscutoff_rna, padj < padj.cutoff)

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
rownames(heat_annot.rna) <- colnames(all50_norm.rna)
heat_annot.rpf <- sampleTable.rpf %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select("Cell Line" = cell_line, "Treatment" = treatment)

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
# heat_rownames_pos.rna <- pos50_norm.rna$gene_id
# heat_rownames_neg.rna <- neg50_norm.rna$gene_id
# heat_rownames_pos.rpf <- pos50_norm.rpf$gene_id
# heat_rownames_neg.rpf <- neg50_norm.rpf$gene_id

heat_rownames_all.rna <- all50_norm.rna$symbol
heat_rownames_all.rpf <- all50_norm.rpf$symbol

# Run pheatmap using the metadata data frame for the annotation
# A <- pheatmap(as.matrix(pos50_norm.rna[2:9]), 
#          color = heat_colors, 
#          cluster_rows = T, 
#          show_rownames = T,
#          labels_row = heat_rownames_pos.rna,
#          annotation = heat_annot.rna, 
#          border_color = NA, 
#          fontsize = 10, 
#          scale = "row", 
#          fontsize_row = 6, 
#          height = 20)
# 
# B <- pheatmap(neg50_norm.rna[2:9], 
#          color = heat_colors, 
#          cluster_rows = T, 
#          show_rownames = T,
#          labels_row = heat_rownames_neg.rna,
#          annotation = heat_annot.rna, 
#          border_color = NA, 
#          fontsize = 10, 
#          scale = "row", 
#          fontsize_row = 6, 
#          height = 20)
# 
# C <- pheatmap(pos50_norm.rpf[2:9], 
#          color = heat_colors, 
#          cluster_rows = T, 
#          show_rownames = T,
#          labels_row = heat_rownames_pos.rpf,
#          annotation = heat_annot.rpf, 
#          border_color = NA, 
#          fontsize = 10, 
#          scale = "row", 
#          fontsize_row = 6, 
#          height = 20)
# 
# D <- pheatmap(neg50_norm.rpf[2:9], 
#          color = heat_colors, 
#          cluster_rows = T, 
#          show_rownames = T,
#          labels_row = heat_rownames_neg.rpf,
#          annotation = heat_annot.rpf, 
#          border_color = NA, 
#          fontsize = 10, 
#          scale = "row", 
#          fontsize_row = 6, 
#          height = 20)

rownames(heat_annot.rna) <- colnames(all50_norm.rna[2:9])
samp_annot.rna <- heat_annot.rna %>% arrange(Treatment)
samp_ordered.rna <- all50_norm.rna[rownames(samp_annot.rna)]

E <- pheatmap(samp_ordered.rna, # Width: 600, Height: 1000
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

rownames(heat_annot.rpf) <- colnames(all50_norm.rpf[2:9])
samp_annot.rpf <- heat_annot.rpf %>% arrange(Treatment)
samp_ordered.rpf <- all50_norm.rpf[rownames(samp_annot.rpf)]

G <- pheatmap(samp_ordered.rpf, 
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
  E[[4]],
  G[[4]],
  align = 'vh',
  labels = "AUTO",
  nrow = 1
)

# Add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(prow)
