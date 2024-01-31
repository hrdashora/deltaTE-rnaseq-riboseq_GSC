### STEP 2: Quality assessment and visualization ###

# Load .RData file from Step 1
## Move .RData file from Isilon to working directory with the following command in the terminal:
# rsync -r -progress -size-only /mnt/isilon/w_stemcell/yuj2lab/HRD/Sequencing/'RNA Core' /Users/dashorh/rnaseq-riboseq/data/

load("data/RNA-RIBO_SummarizedExperiment.RData")

# Load libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(PoiClaClu)
  library(cowplot)
  library(gridExtra)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(scuttle)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(viridis)
})

# Construct a DESeqDataSet object -----------------------------------------
# Confirm that sample information is present
colData(se.rnaseq) 
colData(se.riboseq)

# Determine the design formula from the factors in the sample table
se.rnaseq$treatment
se.riboseq$treatment
se.rnaseq$cell_line
se.riboseq$cell_line

# Set the appropriate reference levels
se.rnaseq$treatment <- relevel(se.rnaseq$treatment, ref = "normoxia")
se.riboseq$treatment <- relevel(se.riboseq$treatment, ref = "normoxia")

# Check the millions of fragments that uniquely align
## Samples from the CLP source and JSY source are duplicated. Indexes of
## duplicated entries in both SummarizedExperiment objects are shown here.
round(colSums(assay(se.rnaseq)) / 1e6, 1)
round(colSums(assay(se.rnaseq)) / 1e6, 1)[1:4] ==
  round(colSums(assay(se.rnaseq)) / 1e6, 1)[c(16:17,20:21)]

round(colSums(assay(se.riboseq)) / 1e6, 1)[1:4] ==
  round(colSums(assay(se.riboseq)) / 1e6, 1)[c(11:12,15:16)]

# Perform two-dimensional subsetting of the SummarizedExperiment objects while
# creating the DESeqDataSet objects
dds.rna_hn <- DESeqDataSet(se.rnaseq[,se.rnaseq$source == "JSY"],
                           design = ~ cell_line + treatment)
dds.rpf_hn <- DESeqDataSet(se.riboseq[,se.riboseq$source == "JSY"],
                           design = ~ cell_line + treatment)

dds.rna_gscnpc <- DESeqDataSet(se.rnaseq[,se.rnaseq$treatment == "normoxia"],
                               design = ~ cell_line)
dds.rpf_gscnpc <- DESeqDataSet(se.riboseq[,se.riboseq$treatment == "normoxia"],
                               design = ~ cell_line)

# Apply minimal pre-filtering of the data-set
## Keep only rows that have a count of at least 10 for a minimal number of
## samples in a group
source("prefilter_genes.R")
dds.rna_hn <- prefilter_genes(dds = dds.rna_hn, group = "treatment")
dds.rpf_hn <- prefilter_genes(dds = dds.rpf_hn, group = "treatment")
dds.rna_gscnpc <- prefilter_genes(dds = dds.rna_gscnpc, group = "cell_line")
dds.rpf_gscnpc <- prefilter_genes(dds = dds.rpf_gscnpc, group = "cell_line")

# Visualize the data following count normalization-----------------------------------------------------
# Perform rlog and vst transformations
## RNA-seq with hypoxia versus normoxia contrast
rld.rna_hn <- rlog(dds.rna_hn, blind = FALSE)
vsd.rna_hn <- vst(dds.rna_hn, blind = FALSE)
dds.rna_hn <- estimateSizeFactors(dds.rna_hn) # account for sequencing depth
df.rna_hn <- bind_rows(
  as_tibble(log2(counts(dds.rna_hn, normalized = TRUE)[, 1:2] + 1)) %>% 
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld.rna_hn)[, 1:2]) %>%
    mutate(transformation = "rlog") %>%
    rename_with(~ c('V1', 'V2'), 1:2),
  as_tibble(assay(vsd.rna_hn)[, 1:2]) %>%
    mutate(transformation = "vst") %>%
    rename_with(~ c('V1', 'V2'), 1:2))
colnames(df.rna_hn)[1:2] <- c("x", "y")
ggplot(df.rna_hn, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid(. ~ transformation)
## RIBO-seq with hypoxia versus normoxia contrast
rld.rpf_hn <- rlog(dds.rpf_hn, blind = FALSE)
vsd.rpf_hn <- vst(dds.rpf_hn, blind = FALSE)
dds.rpf_hn <- estimateSizeFactors(dds.rpf_hn) # account for sequencing depth
df.rpf_hn <- bind_rows(
  as_tibble(log2(counts(dds.rpf_hn, normalized = TRUE)[, 1:2] + 1)) %>% 
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld.rpf_hn)[, 1:2]) %>%
    mutate(transformation = "rlog") %>%
    rename_with(~ c('V1', 'V2'), 1:2),
  as_tibble(assay(vsd.rpf_hn)[, 1:2]) %>%
    mutate(transformation = "vst") %>%
    rename_with(~ c('V1', 'V2'), 1:2))
colnames(df.rpf_hn)[1:2] <- c("x", "y")
ggplot(df.rpf_hn, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid(. ~ transformation)
## RNA-seq with gsc versus npc contrast
rld.rna_gscnpc <- rlog(dds.rna_gscnpc, blind = FALSE)
vsd.rna_gscnpc <- vst(dds.rna_gscnpc, blind = FALSE)
dds.rna_gscnpc <- estimateSizeFactors(dds.rna_gscnpc) # account for sequencing depth
## RIBO-seq with gsc versus npc contrast
rld.rpf_gscnpc <- rlog(dds.rpf_gscnpc, blind = FALSE)
vsd.rpf_gscnpc <- vst(dds.rpf_gscnpc, blind = FALSE)
dds.rpf_gscnpc <- estimateSizeFactors(dds.rpf_gscnpc) # account for sequencing depth

# Calculate sample distances to measure similarities and dissimilarities between samples
rna_hn_sampleDists <- dist(t(assay(vsd.rna_hn)))
rpf_hn_sampleDists <- dist(t(assay(vsd.rpf_hn)))
rna_gscnpc_sampleDists <- dist(t(assay(vsd.rna_gscnpc)))
rpf_gscnpc_sampleDists <- dist(t(assay(vsd.rpf_gscnpc)))

# Plot heatmap of sample distances using vst-transformed values
rna_hn_sampleDistMatrix <- as.matrix(rna_hn_sampleDists)
rownames(rna_hn_sampleDistMatrix) <- paste(vsd.rna_hn$cell_line,
                                           vsd.rna_hn$treatment,
                                           sep = "-")
colnames(rna_hn_sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(rna_hn_sampleDistMatrix, 
         clustering_distance_rows = rna_hn_sampleDists,
         clustering_distance_cols = rna_hn_sampleDists,
         col = colors)

rpf_hn_sampleDistMatrix <- as.matrix(rpf_hn_sampleDists)
rownames(rpf_hn_sampleDistMatrix) <- paste(vsd.rpf_hn$cell_line,
                                           vsd.rpf_hn$treatment,
                                           sep = "-")
colnames(rpf_hn_sampleDistMatrix) <- NULL
pheatmap(rpf_hn_sampleDistMatrix, 
         clustering_distance_rows = rpf_hn_sampleDists,
         clustering_distance_cols = rpf_hn_sampleDists,
         col = colors)

rna_gscnpc_sampleDistMatrix <- as.matrix(rna_gscnpc_sampleDists)
rownames(rna_gscnpc_sampleDistMatrix) <- paste(vsd.rna_gscnpc$cell_line,
                                               vsd.rna_gscnpc$replicate,
                                               sep = "-")
colnames(rna_gscnpc_sampleDistMatrix) <- NULL
pheatmap(rna_gscnpc_sampleDistMatrix, 
         clustering_distance_rows = rna_gscnpc_sampleDists,
         clustering_distance_cols = rna_gscnpc_sampleDists,
         col = colors)
# Remove 318 Replicate 1 and 2 from further analysis due to inter-sample dissimilarity
dds.rna_gscnpc <- dds.rna_gscnpc[,3:15]
# Re-transform counts from the DESeqDataSet object, sample distances are not
# utilized further so updates not needed
rld.rna_gscnpc <- rlog(dds.rna_gscnpc, blind = FALSE)
vsd.rna_gscnpc <- vst(dds.rna_gscnpc, blind = FALSE)
dds.rna_gscnpc <- estimateSizeFactors(dds.rna_gscnpc)

rpf_gscnpc_sampleDistMatrix <- as.matrix(rpf_gscnpc_sampleDists)
rownames(rpf_gscnpc_sampleDistMatrix) <- paste(vsd.rpf_gscnpc$cell_line,
                                               vsd.rpf_gscnpc$replicate,
                                               sep = "-")
colnames(rpf_gscnpc_sampleDistMatrix) <- NULL
pheatmap(rpf_gscnpc_sampleDistMatrix, 
         clustering_distance_rows = rpf_gscnpc_sampleDists,
         clustering_distance_cols = rpf_gscnpc_sampleDists,
         col = colors) 


# Generate PCA plots
## RNA-seq with hypoxia versus normoxia contrast
rna_hn_pcaData <- plotPCA(vsd.rna_hn, intgroup = c("treatment", "cell_line"), returnData = TRUE)
rna_hn_percentVar <- round(100 * attr(rna_hn_pcaData, "percentVar"))
rna_hn_pcaPlot <- ggplot(rna_hn_pcaData, aes(x = PC1, y = PC2, color = treatment, shape = cell_line)) +
  geom_point(size = 3) +
  labs(color = "Treatment", shape = "Cell Line") +
  xlab(paste0("PC1: ",rna_hn_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",rna_hn_percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
rna_hn_pcaPlot # Publish using Width: 600, Height: 400

## RIBO-seq with hypoxia versus normoxia contrast
rpf_hn_pcaData <- plotPCA(vsd.rpf_hn, intgroup = c("treatment", "cell_line"), returnData = TRUE)
rpf_hn_percentVar <- round(100 * attr(rpf_hn_pcaData, "percentVar"))
rpf_hn_pcaPlot <- ggplot(rpf_hn_pcaData, aes(x = PC1, y = PC2, color = treatment, shape = cell_line)) +
  geom_point(size = 3) +
  labs(color = "Treatment", shape = "Cell Line") +
  xlab(paste0("PC1: ",rpf_hn_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",rpf_hn_percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
rpf_hn_pcaPlot # Publish using Width: 600, Height: 400

legend <- cowplot::get_legend( # extract the legend from one of the plots
  rna_hn_pcaPlot + theme(legend.box.margin = margin(0, 0, 0, 12))) 
prow <- plot_grid( # arrange the plots in a single row
  rna_hn_pcaPlot + theme(legend.position = "none"),
  rpf_hn_pcaPlot + theme(legend.position = "none"),
  align = 'vh',
  nrow = 1)
fig <- plot_grid(prow, legend, rel_widths = c(3, .4)) # add the legend to the row made earlier
grid.arrange(fig)

## RNA-seq with gsc versus npc contrast
rna_gscnpc_pcaData <- plotPCA(vsd.rna_gscnpc, intgroup = "cell_line", returnData = TRUE)
rna_gscnpc_percentVar <- round(100 * attr(rna_gscnpc_pcaData, "percentVar"))
rna_gscnpc_pcaPlot <- ggplot(rna_gscnpc_pcaData, aes(x = PC1, y = PC2, color = cell_line)) +
  geom_point(size = 3) +
  labs(color = "Cell Line") +
  xlab(paste0("PC1: ",rna_gscnpc_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",rna_gscnpc_percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
rna_gscnpc_pcaPlot # Publish using Width: 600, Height: 400

## RIBO-seq with gsc versus npc contrast
rpf_gscnpc_pcaData <- plotPCA(vsd.rpf_gscnpc, intgroup = "cell_line", returnData = TRUE)
rpf_gscnpc_percentVar <- round(100 * attr(rpf_gscnpc_pcaData, "percentVar"))
rpf_gscnpc_pcaPlot <- ggplot(rpf_gscnpc_pcaData, aes(x = PC1, y = PC2, color = cell_line)) +
  geom_point(size = 3) +
  labs(color = "Cell Line") +
  xlab(paste0("PC1: ",rpf_gscnpc_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",rpf_gscnpc_percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
rpf_gscnpc_pcaPlot # Publish using Width: 600, Height: 400

legend <- cowplot::get_legend( # extract the legend from one of the plots
  rna_gscnpc_pcaPlot + theme(legend.box.margin = margin(0, 0, 0, 12)))
prow <- plot_grid( # arrange the plots in a single row
  rna_gscnpc_pcaPlot + theme(legend.position = "none"),
  rpf_gscnpc_pcaPlot + theme(legend.position = "none"),
  align = 'vh',
  nrow = 1)
fig <- plot_grid(prow, legend, rel_widths = c(3, .4)) # add the legend to the row made earlier
grid.arrange(fig)
