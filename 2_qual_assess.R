### STEP 2: Quality assessment and visualization ###
load("RNA-RIBO_SummarizedExperiment.RData")

# Load libraries
library(DESeq2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(vsn)
library(cowplot)
library(gridExtra)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(scuttle)
library(tidyr)
library(tibble)
library(stringr)
library(viridis)

# library(SummarizedExperiment)
# library(MASS)
# library(GenomicFeatures)
# library(GenomicRanges)

# Determine the design formula
se.rna$treatment
se.rpf$treatment

se.rna$treatment <- relevel(se.rna$treatment, ref = "Normoxic")
se.rpf$treatment <- relevel(se.rpf$treatment, ref = "Normoxic")

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

# countdata.rna <- fc.rna$counts # starting from count matrices
# head(countdata.rna, 3)
# coldata.rna <- sampleTable.rna
# 
# ddsMat.rna <- DESeqDataSetFromMatrix(countData = countdata.rna,
#                                      colData = coldata.rna,
#                                      design = ~ cell_line + treatment)
# 
# countdata.rpf <- fc.rpf$counts
# head(countdata.rpf, 3)
# coldata.rpf <- sampleTable.rpf
# 
# ddsMat.rpf <- DESeqDataSetFromMatrix(countData = countdata.rpf,
#                                      colData = coldata.rpf,
#                                      design = ~ cell_line + treatment)

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
rna_pcaPlot # Width: 600, Height: 400

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

Figure1 <- plot_grid(prow, legend, rel_widths = c(3, .4)) # add the legend to the row made earlier
#png("~/Desktop/Figure1.png", width = 20, height = 10, units = 'cm', res = 300)
grid.arrange(Figure1)

# MDS plot
rna_mds <- as.data.frame(colData(rld.rna)) %>%
  cbind(cmdscale(rna_sampleDistMatrix))
ggplot(rna_mds, aes(x = `1`, y = `2`, color = treatment, shape = cell_line)) +
  geom_point(size = 3) + coord_fixed()

rna_mdsPois <- as.data.frame(colData(dds.rna)) %>%
  cbind(cmdscale(rna_samplePoisDistMatrix))
ggplot(rna_mdsPois, aes(x = `1`, y = `2`, color = treatment, shape = cell_line)) +
  geom_point(size = 3) + coord_fixed()

rpf_mds <- as.data.frame(colData(rld.rpf)) %>%
  cbind(cmdscale(rpf_sampleDistMatrix))
ggplot(rpf_mds, aes(x = `1`, y = `2`, color = treatment, shape = cell_line)) +
  geom_point(size = 3) + coord_fixed()

rpf_mdsPois <- as.data.frame(colData(dds.rpf)) %>%
  cbind(cmdscale(rpf_samplePoisDistMatrix))
ggplot(rpf_mdsPois, aes(x = `1`, y = `2`, color = treatment, shape = cell_line)) +
  geom_point(size = 3) + coord_fixed()

# Calculate TPMs
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by = "gene")
exonicLengths <- sum(width(GenomicRanges::reduce(ebg)))
rna_keep <- intersect(names(exonicLengths), rownames(dds.rna))
rpf_keep <- intersect(names(exonicLengths), rownames(dds.rpf))

rna_tpm <- calculateTPM(dds.rna[rna_keep,], lengths = exonicLengths[rna_keep])
rpf_tpm <- calculateTPM(dds.rpf[rpf_keep,], lengths = exonicLengths[rpf_keep])

# # Alternative TPM Calculation
# exonicLengths_kilobases <- exonicLengths/1000 # convert to kilobases
# rpk <- sweep(assay(se.rna[keep,]),1,exonicLengths_kilobases[keep],'/')
# pmScaling <- colSums(rpk)/1e6 # "per million" scaling factor
# tpm_manual <- sweep(rpk, 2, pmScaling, '/')

# Create scatterplot of TPM values
rna_tpm <- rna_tpm %>%
  as_tibble(rownames = "gene_id") %>% # convert matrix to dataframe
  rename_with(~ sampleTable.rna$sample, starts_with('V')) %>% # adjust column names in tpm dataframe
  pivot_longer(cols = 2:9,
               names_to = "sample",
               values_to = "tpm") %>% # convert from wide to long data format
  left_join(sampleTable.rna, by = "sample")

rpf_tpm <- rpf_tpm %>%
  as_tibble(rownames = "gene_id") %>% # convert matrix to dataframe
  rename_with(~ sampleTable.rpf$sample, starts_with('V')) %>% # adjust column names in tpm dataframe
  pivot_longer(cols = 2:9,
               names_to = "sample",
               values_to = "tpm") %>% # convert from wide to long data format
  left_join(sampleTable.rpf, by = "sample") %>%
  mutate(sample = str_remove(sample, 'r'))

all_tpm <- merge(rna_tpm,
                 rpf_tpm,
                 by = c("gene_id", "sample", "cell_line", "treatment"),
                 suffixes = c("_rna","_rpf")) %>% 
  dplyr::select(-sample)
  
byCondition_tpm <- all_tpm %>% group_by(gene_id, treatment) %>% # calculate average TPM across replicate and cell line
  summarise(mean_rna = mean(tpm_rna),
            mean_rpf = mean(tpm_rpf)) %>%
  mutate(log2.rna = log(mean_rna + 1, 2),
         log2.rpf = log(mean_rpf + 1, 2))

mod <- lm(log2.rpf ~ log2.rna, data = byCondition_tpm)
ggplot(byCondition_tpm, aes(log2.rna, log2.rpf)) +
  geom_point(aes(color = treatment),
             alpha = 1/10) +
  geom_abline(intercept = coef(mod)[1], slope = coef(mod)[2], col = 'black', linetype = 3) +
  labs(x = "RNA-seq (log2 TPM)",
       y = "RIBO-seq (log2 TPM)",
       color = "Treatment") +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 3.5, vjust = 2, parse = TRUE,
           label = paste0("R^2==", round(summary(mod)$r.squared,3))) +
  theme_classic()

# Create density plot
byGene_tpm <- dplyr::select(byCondition_tpm,-log2.rna,-log2.rpf) %>%
  pivot_wider(names_from = treatment, values_from = c(mean_rna, mean_rpf)) %>%
  mutate(hypo_log2_rna = log(mean_rna_Hypoxic + 1, 2),
         norm_log2_rna = log(mean_rna_Normoxic + 1, 2),
         hypo_log2_rpf = log(mean_rpf_Hypoxic + 1, 2),
         norm_log2_rpf = log(mean_rpf_Normoxic + 1, 2))

source("get_density.R")
set.seed(2022)

byGene_tpm$density_rna <- get_density(byGene_tpm$norm_log2_rna, byGene_tpm$hypo_log2_rna, h = c(3,3), n = 100)
byGene_tpm$density_rpf <- get_density(byGene_tpm$norm_log2_rpf, byGene_tpm$hypo_log2_rpf, h = c(3,3), n = 100)

mod_rna <- lm(hypo_log2_rna ~ norm_log2_rna, data = byGene_tpm)
tpmPlot_rna <- ggplot(byGene_tpm) + geom_point(aes(norm_log2_rna, hypo_log2_rna, color = density_rna)) +
  scale_color_viridis() +
  geom_abline(intercept = coef(mod_rna)[1], slope = coef(mod_rna)[2], col='black', linetype = 3) +
  labs(x = "Normoxia (log2 TPM)",
       y = "Hypoxia (log2 TPM)",
       color = "Density") +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 3.5, vjust = 2, parse=TRUE,
           label = paste0("R^2==", round(summary(mod_rna)$r.squared,3))) +
  theme_classic() +
  theme(legend.justification = c(1,0), legend.position = c(0.95,0.05), legend.key.size = unit(3, 'mm'))
tpmPlot_rna

mod_rpf <- lm(hypo_log2_rpf ~ norm_log2_rpf, data = byGene_tpm)
tpmPlot_rpf <- ggplot(byGene_tpm) + geom_point(aes(norm_log2_rpf, hypo_log2_rpf, color = density_rpf)) +
  scale_color_viridis() +
  geom_abline(intercept = coef(mod_rpf)[1], slope = coef(mod_rpf)[2], col='black', linetype = 3) +
  labs(x = "Normoxia (log2 TPM)",
       y = "Hypoxia (log2 TPM)",
       color = "Density") +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 3.5, vjust = 2, parse=TRUE,
           label = paste0("R^2==", round(summary(mod_rpf)$r.squared,3))) +
  theme_classic() +
  theme(legend.justification = c(1,0), legend.position = c(0.95,0.05), legend.key.size = unit(3, 'mm'))
tpmPlot_rpf

# Combine QC plots into multi-panel figure
plot_grid(
  rna_pcaPlot,
  rpf_pcaPlot,
  tpmPlot_rna,
  tpmPlot_rpf,
  align = 'h',
  nrow = 2,
  ncol = 2,
  labels = "AUTO",
  hjust = -1
)

# Calculate translational efficiency
te <- byCondition_tpm %>%
  mutate(log2_te = log((mean_rpf) + 1, 2) - log((mean_rna) + 1, 2)) %>%
  mutate(te = 2^log2_te) %>%
  group_by(gene_id)

delta_te <- summarise(te,
                      log2_delta_te = log((mean_rpf[treatment == "Hypoxic"]) + 1, 2) # TE_Hypoxia/TE_Normoxia
                      - log((mean_rna[treatment == "Hypoxic"]) + 1, 2)
                      - log((mean_rpf[treatment == "Normoxic"]) + 1, 2)
                      + log((mean_rna[treatment == "Normoxic"]) + 1, 2)) %>%
  mutate(delta_te = 2^log2_delta_te)

ggplot(delta_te, aes(log2_delta_te,after_stat(density))) + geom_freqpoly(binwidth = 0.1)

## MOLLY
library(ggpubr)
library(org.Hs.eg.db)
str_rna <- rna_tpm$gene_id
rna_tpm$symbol <- mapIds(org.Hs.eg.db,
                         keys = str_rna,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")

# plot_data <- rna_tpm %>%
#   group_by(cell_line, symbol) %>%
#   filter(treatment == "Normoxic",
#          symbol %in% c("SMARCA1","SMARCA5","CHD4", "ACTB")) %>%
#   summarise(avg_tpm = signif(log2(mean(tpm)+1), digits = 3))

plot_data <- rna_tpm %>%
  group_by(cell_line, symbol) %>%
  filter(treatment == "Normoxic",
         symbol %in% c("LILRB4","CD276","ROR1")) %>%
  summarise(avg_tpm = signif(log2(mean(tpm)+1), digits = 3))

ggbarplot(plot_data, "symbol", "avg_tpm",
          fill = "cell_line", color = "cell_line", palette = "Paired",
          label = TRUE,
          position = position_dodge(0.8),
          xlab = "Gene Symbol",
          ylab = "log2(TPM + 1)") +
  labs(fill = "GSC Line", color = "GSC Line") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

