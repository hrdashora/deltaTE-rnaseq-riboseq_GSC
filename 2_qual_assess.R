### STEP 2: Quality assessment and visualization ###

## This script will

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

# Determine the design formula
se.rna$treatment
se.rpf$treatment

se.rna$treatment <- relevel(se.rna$treatment, ref = "Normoxia")
se.rpf$treatment <- relevel(se.rpf$treatment, ref = "Normoxia")

se.rna$cell_line
se.rpf$cell_line

# Construct a DESeqDataSet object
round(colSums(assay(se.rna)) / 1e6, 1) # check the millions of fragments that uniquely align
round(colSums(assay(se.rpf)) / 1e6, 1)
colData(se.rna) # confirm that sample information is present
colData(se.rpf)
dds.rna <- DESeqDataSet(se.rna, design = ~ cell_line + treatment) # starting from SummarizedExperiment
dds.rpf <- DESeqDataSet(se.rpf, design = ~ cell_line + treatment)

# Apply minimal pre-filtering of the data-set
nrow(dds.rna)
nrow(dds.rpf)

dds.rna <- dds.rna[ rowSums(counts(dds.rna)) > 1, ]
dds.rpf <- dds.rpf[ rowSums(counts(dds.rpf)) > 1, ]

nrow(dds.rna)
nrow(dds.rpf)

# Rlog and vst transformations
rld.rna <- rlog(dds.rna, blind = FALSE)
head(assay(rld.rna), 3)
meanSdPlot(assay(rld.rna), ranks = FALSE)
vsd.rna <- vst(dds.rna, blind = FALSE)
head(assay(vsd.rna), 3)
dds.rna <- estimateSizeFactors(dds.rna) # account for sequencing depth

df.rna <- bind_rows(
  as_tibble(log2(counts(dds.rna, normalized = TRUE)[, 1:2] + 1)) %>% 
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld.rna)[, 1:2]) %>%
    mutate(transformation = "rlog") %>%
    rename_with(~ c('V1', 'V2'), 1:2),
  as_tibble(assay(vsd.rna)[, 1:2]) %>%
    mutate(transformation = "vst") %>%
    rename_with(~ c('V1', 'V2'), 1:2))

colnames(df.rna)[1:2] <- c("x", "y")
ggplot(df.rna, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid(. ~ transformation)

rld.rpf <- rlog(dds.rpf, blind = FALSE)
head(assay(rld.rpf), 3)
meanSdPlot(assay(rld.rpf), ranks = FALSE)
vsd.rpf <- vst(dds.rpf, blind = FALSE)
head(assay(vsd.rpf), 3)
dds.rpf <- estimateSizeFactors(dds.rpf) # account for sequencing depth

df.rpf <- bind_rows(
  as_tibble(log2(counts(dds.rpf, normalized = TRUE)[, 1:2] + 1)) %>% 
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld.rpf)[, 1:2]) %>%
    mutate(transformation = "rlog") %>%
    rename_with(~ c('V1', 'V2'), 1:2),
  as_tibble(assay(vsd.rpf)[, 1:2]) %>%
    mutate(transformation = "vst") %>%
    rename_with(~ c('V1', 'V2'), 1:2))

colnames(df.rpf)[1:2] <- c("x", "y")
ggplot(df.rpf, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid(. ~ transformation)

# Calculate sample distances
rna_sampleDists <- dist(t(assay(rld.rna)))
rpf_sampleDists <- dist(t(assay(rld.rpf)))

rna_sampleDistMatrix <- as.matrix(rna_sampleDists)
rownames(rna_sampleDistMatrix) <- paste(rld.rna$sample)
colnames(rna_sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(rna_sampleDistMatrix, # heatmap of sample distances using rlog-transformed values
         clustering_distance_rows = rna_sampleDists,
         clustering_distance_cols = rna_sampleDists,
         col = colors)

rpf_sampleDistMatrix <- as.matrix(rpf_sampleDists)
rownames(rpf_sampleDistMatrix) <- paste(rld.rpf$sample)
colnames(rpf_sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(rpf_sampleDistMatrix, # heatmap of sample distances using rlog-transformed values
         clustering_distance_rows = rpf_sampleDists,
         clustering_distance_cols = rpf_sampleDists,
         col = colors)

rna_poisd <- PoissonDistance(t(counts(dds.rna)))
rna_samplePoisDistMatrix <- as.matrix(rna_poisd$dd)
rownames(rna_samplePoisDistMatrix) <- paste(rld.rna$sample)
colnames(rna_samplePoisDistMatrix) <- NULL
pheatmap(rna_samplePoisDistMatrix, # heatmap of the sample distances using the Poisson distance
         clustering_distance_rows = rna_poisd$dd,
         clustering_distance_cols = rna_poisd$dd,
         col = colors)

rpf_poisd <- PoissonDistance(t(counts(dds.rpf)))
rpf_samplePoisDistMatrix <- as.matrix(rpf_poisd$dd)
rownames(rpf_samplePoisDistMatrix) <- paste(rld.rpf$sample)
colnames(rpf_samplePoisDistMatrix) <- NULL
pheatmap(rpf_samplePoisDistMatrix, # heatmap of the sample distances using the Poisson distance
         clustering_distance_rows = rpf_poisd$dd,
         clustering_distance_cols = rpf_poisd$dd,
         col = colors)

# PCA plot
rna_pcaData <- plotPCA(rld.rna, intgroup = c("treatment", "cell_line"), returnData = TRUE)
rna_percentVar <- round(100 * attr(rna_pcaData, "percentVar"))
rna_pcaPlot <- ggplot(rna_pcaData, aes(x = PC1, y = PC2, color = treatment, shape = cell_line)) +
  geom_point(size = 3) +
  labs(color = "Treatment", shape = "Cell Line") +
  xlab(paste0("PC1: ",rna_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",rna_percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
rna_pcaPlot # Publish using Width: 600, Height: 400

rpf_pcaData <- plotPCA(rld.rpf, intgroup = c("treatment", "cell_line"), returnData = TRUE)
rpf_percentVar <- round(100 * attr(rpf_pcaData, "percentVar"))
rpf_pcaPlot <- ggplot(rpf_pcaData, aes(x = PC1, y = PC2, color = treatment, shape = cell_line)) +
  geom_point(size = 3) +
  labs(color = "Treatment", shape = "Cell Line") +
  xlab(paste0("PC1: ",rpf_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",rpf_percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
rpf_pcaPlot

legend <- cowplot::get_legend( # extract the legend from one of the plots
  rna_pcaPlot + theme(legend.box.margin = margin(0, 0, 0, 12))) 

prow <- plot_grid( # arrange the plots in a single row
  rna_pcaPlot + theme(legend.position = "none"),
  rpf_pcaPlot + theme(legend.position = "none"),
  align = 'vh',
  nrow = 1)

fig <- plot_grid(prow, legend, rel_widths = c(3, .4)) # add the legend to the row made earlier
grid.arrange(fig)
