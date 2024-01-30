### STEP 3: Calculate TPM and visualize TE (translational efficiency) ###

# Calculate TPMs
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by = "gene")
exonicLengths <- sum(width(GenomicRanges::reduce(ebg)))
rna_keep <- intersect(names(exonicLengths), rownames(dds.rna))
rpf_keep <- intersect(names(exonicLengths), rownames(dds.rpf))

rna_tpm <- calculateTPM(dds.rna[rna_keep,], lengths = exonicLengths[rna_keep]) # returns a matrix: N (21,161) x M (8)
rpf_tpm <- calculateTPM(dds.rpf[rpf_keep,], lengths = exonicLengths[rpf_keep])

# Create scatterplot of TPM values
rna_tpm <- rna_tpm %>%
  as_tibble(rownames = "gene_id") %>% # convert matrix to dataframe
  rename_with(~ sampleTable.rna$sample, .cols = starts_with('V')) %>% # adjust column names in tpm dataframe
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

byCondition_tpm <- all_tpm %>%
  group_by(gene_id, treatment) %>% # calculate average TPM across replicates and cell line
  summarise(mean_rna = mean(tpm_rna),
            mean_rpf = mean(tpm_rpf)) %>%
  mutate(log2_rna = log(mean_rna + 1, 2),
         log2_rpf = log(mean_rpf + 1, 2))

mod <- lm(log2_rpf ~ log2_rna, data = byCondition_tpm)
ggplot(byCondition_tpm, aes(log2_rna, log2_rpf)) +
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
byGene_tpm <- dplyr::select(byCondition_tpm,-log2_rna,-log2_rpf) %>%
  pivot_wider(names_from = treatment, values_from = c(mean_rna, mean_rpf)) %>%
  mutate(hypo_log2_rna = log(mean_rna_Hypoxia + 1, 2),
         norm_log2_rna = log(mean_rna_Normoxia + 1, 2),
         hypo_log2_rpf = log(mean_rpf_Hypoxia + 1, 2),
         norm_log2_rpf = log(mean_rpf_Normoxia + 1, 2))

source("get_density.R")
set.seed(2023)

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
## NOTE: When calculating the Log2 transform of TPM values, use Log2(TPM + 1).
## and only apply that calculation to the base TPM values, not to secondary ratios
te <- byCondition_tpm %>%
  mutate(log2_te = log((mean_rpf) + 1, 2) - log((mean_rna) + 1, 2)) %>%
  mutate(te = 2^log2_te) %>%
  group_by(gene_id)

delta_te <- summarise(te,
                      log2_delta_te = log2_te[treatment == "Hypoxia"] 
                      - log2_te[treatment == "Normoxia"]) %>%
  mutate(delta_te = 2^log2_delta_te) # TE_Hypoxia/TE_Normoxia

ggplot(delta_te, aes(log2_delta_te, after_stat(density))) + geom_freqpoly(binwidth = 0.1)
