### STEP 3: Calculate TPM and visualize TE (translational efficiency) ###

## This script uses the DESeq2 objects prepared in Step 2 to calculate
## TPM (transcripts per million) and visualize TE (translational efficiency) 
## using scatterplots and density plots.

# Run the script in Step 2 in its entirely to generate the necessary variables

# Calculate TPMs ----------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by = "gene")
exonicLengths <- sum(width(GenomicRanges::reduce(ebg)))

rna_hn_keep <- intersect(names(exonicLengths), rownames(dds.rna_hn))
rpf_hn_keep <- intersect(names(exonicLengths), rownames(dds.rpf_hn))
rna_gscnpc_keep <- intersect(names(exonicLengths), rownames(dds.rna_gscnpc))
rpf_gscnpc_keep <- intersect(names(exonicLengths), rownames(dds.rpf_gscnpc))

rna_hn_tpm <- scuttle::calculateTPM(dds.rna_hn[rna_hn_keep,],
                                    lengths = exonicLengths[rna_hn_keep]) # returns a matrix
rpf_hn_tpm <- scuttle::calculateTPM(dds.rpf_hn[rpf_hn_keep,],
                                    lengths = exonicLengths[rpf_hn_keep])
rna_gscnpc_tpm <- scuttle::calculateTPM(dds.rna_gscnpc[rna_gscnpc_keep,],
                                       lengths = exonicLengths[rna_gscnpc_keep])
rpf_gscnpc_tpm <- scuttle::calculateTPM(dds.rpf_gscnpc[rpf_gscnpc_keep,],
                                       lengths = exonicLengths[rpf_gscnpc_keep])

# Create merged data frame of TPM values
rna_hn_tpm <- rna_hn_tpm %>%
  as_tibble(rownames = "gene_id") %>% # convert matrix to dataframe
  rename_with(~ colData(dds.rna_hn)$file_name, .cols = starts_with('V')) %>% # adjust column names in tpm dataframe
  pivot_longer(cols = 2:9,
               names_to = "file_name",
               values_to = "tpm") %>% # convert from wide to long data format
  left_join(colData(dds.rna_hn), by = "file_name", copy = TRUE)

rpf_hn_tpm <- rpf_hn_tpm %>%
  as_tibble(rownames = "gene_id") %>% # convert matrix to dataframe
  rename_with(~ colData(dds.rpf_hn)$file_name, .cols = starts_with('V')) %>% # adjust column names in tpm dataframe
  pivot_longer(cols = 2:9,
               names_to = "file_name",
               values_to = "tpm") %>% # convert from wide to long data format
  left_join(colData(dds.rpf_hn), by = "file_name", copy = TRUE)

all_hn_tpm <- merge(rna_hn_tpm,
                    rpf_hn_tpm,
                    by = c("gene_id", "cell_line", "treatment", "replicate"),
                    suffixes = c("_rna","_rpf")) %>%
  mutate(log2_tpm_rna = log(tpm_rna + 1, 2),
         log2_tpm_rpf = log(tpm_rpf + 1, 2)) %>%
  pivot_wider(id_cols =  c("gene_id", "cell_line", "treatment"),
              names_from = c("replicate"),
              names_prefix = "rep",
              values_from = c("tpm_rna", "tpm_rpf", "log2_tpm_rna", "log2_tpm_rpf")) %>%
  dplyr::select(-starts_with(c("file_name", "source", "sizeFactor"))) # remove unnecessary columns

rna_gscnpc_tpm <- rna_gscnpc_tpm %>%
  as_tibble(rownames = "gene_id") %>% # convert matrix to dataframe
  rename_with(~ colData(dds.rna_gscnpc)$file_name, .cols = starts_with('V')) %>% # adjust column names in tpm dataframe
  pivot_longer(cols = 2:14,
               names_to = "file_name",
               values_to = "tpm") %>% # convert from wide to long data format
  left_join(colData(dds.rna_gscnpc), by = "file_name", copy = TRUE)

rpf_gscnpc_tpm <- rpf_gscnpc_tpm %>%
  as_tibble(rownames = "gene_id") %>% # convert matrix to dataframe
  rename_with(~ colData(dds.rpf_gscnpc)$file_name, .cols = starts_with('V')) %>% # adjust column names in tpm dataframe
  pivot_longer(cols = 2:11,
               names_to = "file_name",
               values_to = "tpm") %>% # convert from wide to long data format
  left_join(colData(dds.rpf_gscnpc), by = "file_name", copy = TRUE)

all_gscnpc_tpm <- merge(rna_gscnpc_tpm, # only replicates present in both sequencing datasets will be merged
                        rpf_gscnpc_tpm,
                        by = c("gene_id", "cell_line", "replicate"),
                        suffixes = c("_rna","_rpf")) %>%
  mutate(log2_tpm_rna = log(tpm_rna + 1, 2),
         log2_tpm_rpf = log(tpm_rpf + 1, 2)) %>%
  pivot_wider(id_cols =  c("gene_id", "cell_line"),
              names_from = c("replicate"),
              names_prefix = "rep",
              values_from = c("tpm_rna", "tpm_rpf", "log2_tpm_rna", "log2_tpm_rpf")) %>%
  dplyr::select(-starts_with(c("file_name", "source", "treatment", "sizeFactor"))) # remove unnecessary columns


# Create facet grid scatter plots of TPMs for each replicate
treat.labs <- c("Normoxia", "Hypoxia") # New facet label names for treatment variable
names(treat.labs) <- c("normoxia", "hypoxia")

ggplot(all_hn_tpm, aes(log2_tpm_rna_rep1, log2_tpm_rna_rep2)) +
  geom_point() +
  facet_grid(cols = vars(cell_line),
             rows = vars(treatment),
             labeller = labeller(treatment = treat.labs)) +
  geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 3) +
  labs(x = "Replicate 1 (log2 TPM)",
       y = "Replicate 2 (log2 TPM)",
       title = "RNA-seq TPM values")

ggplot(all_hn_tpm, aes(log2_tpm_rpf_rep1, log2_tpm_rpf_rep2)) +
  geom_point() +
  facet_grid(cols = vars(cell_line),
             rows = vars(treatment),
             labeller = labeller(treatment = treat.labs)) +
  geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 3) +
  labs(x = "Replicate 1 (log2 TPM)",
       y = "Replicate 2 (log2 TPM)",
       title = "RIBO-seq TPM values")

ggplot(all_gscnpc_tpm, aes(log2_tpm_rna_rep1, log2_tpm_rna_rep2)) +
  geom_point() +
  facet_grid(cols = vars(cell_line)) +
  geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 3) +
  labs(x = "Replicate 1 (log2 TPM)",
       y = "Replicate 2 (log2 TPM)",
       title = "RNA-seq TPM values")

ggplot(all_gscnpc_tpm, aes(log2_tpm_rpf_rep1, log2_tpm_rpf_rep2)) +
  geom_point() +
  facet_grid(cols = vars(cell_line)) +
  geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 3) +
  labs(x = "Replicate 1 (log2 TPM)",
       y = "Replicate 2 (log2 TPM)",
       title = "RIBO-seq TPM values")


# Calculate translational efficiency --------------------------------------
## NOTE: When calculating the Log2 transform of TPM values, use Log2(TPM + 1).
## and only apply that calculation to the base TPM values, not to secondary ratios

delta_te_hn <- all_hn_tpm %>%
  mutate(mean_tpm_rna = rowMeans(dplyr::select(., starts_with("tpm_rna")), na.rm = TRUE),
         mean_tpm_rpf = rowMeans(dplyr::select(., starts_with("tpm_rpf")), na.rm = TRUE)) %>%
  mutate(log2_te = log(mean_tpm_rpf + 1, 2) - log(mean_tpm_rna + 1, 2)) %>%
  mutate(te = 2^log2_te) %>%
  group_by(gene_id, cell_line) %>%
  reframe(log2_delta_te = log2_te[treatment == "hypoxia"] - log2_te[treatment == "normoxia"]) %>%
  mutate(delta_te = 2^log2_delta_te)

ggplot(delta_te_hn, aes(log2_delta_te, after_stat(density))) +
  geom_freqpoly(binwidth = 0.1) +
  facet_grid(cols = vars(cell_line))

te_gscnpc <- all_gscnpc_tpm %>%
  mutate(mean_tpm_rna = rowMeans(dplyr::select(., starts_with("tpm_rna")), na.rm = TRUE),
         mean_tpm_rpf = rowMeans(dplyr::select(., starts_with("tpm_rpf")), na.rm = TRUE)) %>%
  mutate(log2_te = log(mean_tpm_rpf + 1, 2) - log(mean_tpm_rna + 1, 2)) %>%
  mutate(te = 2^log2_te) %>%
  group_by(gene_id, cell_line)
