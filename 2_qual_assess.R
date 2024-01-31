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
dds.rpf_gscnpc <- DESeqDataSet(se.rnaseq[,se.riboseq$treatment == "normoxia"],
                               design = ~ cell_line)

# Apply minimal pre-filtering of the data-set
## Keep only rows that have a count of at least 10 for a minimal number of
## samples
source("prefilter_genes.R")
dds.rna_hn <- prefilter_genes(dds = dds.rna_hn, group = "treatment")


smallestGroupSize <- min(table(dds.rna_hn$treatment))
keep <- rowSums(counts(dds.rna_hn) >= 10) >= smallestGroupSize
dds.rna_hn <- dds.rna_hn[keep,]
nrow(dds.rna_hn)

nrow(dds.rpf_hn)
smallestGroupSize <- min(table(dds.rpf_hn$treatment))
keep <- rowSums(counts(dds.rpf_hn) >= 10) >= smallestGroupSize
dds.rpf_hn <- dds.rpf_hn[keep,]
nrow(dds.rna_hn)

dds.rna_hn$cell_line

nrow(dds.rpf)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

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
