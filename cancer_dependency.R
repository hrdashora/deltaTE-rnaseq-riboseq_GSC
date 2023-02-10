# Load necessary libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(tibble)
library(gridExtra)
library(stringr)
library(depmap)
library(ExperimentHub)

# Load the rnai, crispr, and copyNumber datasets for visualization
## create ExperimentHub query object
eh <- ExperimentHub()
query(eh, "depmap")
rnai <- eh[["EH3080"]]
mutationCalls <- eh[["EH3085"]]
metadata <- eh[["EH3086"]]
TPM <- eh[["EH3084"]]
copyNumber <- eh[["EH3082"]]
crispr <- eh[["EH3081"]]
drug_sensitivity <- eh[["EH3087"]]

## list of dependency scores, highly negative scores imply that a cell line is dependent on that gene
CNS.rnai <- rnai %>% dplyr::select(cell_line, gene_name, dependency) %>%
  dplyr::filter(stringr::str_detect(cell_line, "CENTRAL_NERVOUS")) %>%
  tidyr::drop_na(dependency) %>%
  dplyr::filter(dependency >= -0.50) %>%
  dplyr::arrange(desc(dependency))

CNS.crispr <- crispr %>% dplyr::select(cell_line, gene_name, dependency) %>%
  dplyr::filter(stringr::str_detect(cell_line, "CENTRAL_NERVOUS")) %>%
  tidyr::drop_na(dependency) %>%
  dplyr::filter(dependency >= -0.50) %>%
  dplyr::arrange(desc(dependency))

# Load dependencies enriched in glioma
enrich <- readr::read_csv("data/dependencies enriched in glioma.csv") %>%
  filter(Type == "gene", stringr::str_detect(Dataset, "CRISPR"))

notdep <- table(
  CNS.crispr[which(CNS.crispr$gene_name %in% enrich$`Gene/Compound`),2]) %>%
  as.data.frame() %>%
  arrange(desc(Freq))

te_up <- readr::read_csv("results/te_categories.csv") %>%
  dplyr::filter(category_te == "TE.up") %>% # Add for further filtering: log2_delta_te >= log2(1.75)
  tidyr::drop_na(symbol) %>%
  pull(symbol)

notdep[which(notdep$gene_name %in% te_up),]
targets <- table(CNS.crispr[which(CNS.crispr$gene_name %in% te_up),2]) %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  pull(gene_name)

cyto <- readr::read_delim("results/te_up_reactome_table.txt") %>%
  pull(HitGenes) %>%
  str_split(pattern = ",") %>%
  unlist() %>%
  as_factor()

inter <- intersect(targets, cyto)

# Basic histogram
crispr %>% dplyr::select(gene, gene_name, dependency) %>% 
  dplyr::filter(gene_name == "MCM10") %>% 
  ggplot(aes(x = dependency)) +
  geom_histogram() +
  geom_vline(xintercept = mean(crispr$dependency, na.rm = TRUE),
             linetype = "dotted", color = "red") +
  ggtitle("Histogram of dependency scores for gene MCM10")

meta_crispr <- metadata %>%
  dplyr::select(depmap_id, lineage) %>%
  dplyr::full_join(crispr, by = "depmap_id") %>%
  dplyr::filter(gene_name == "MCM10") %>% 
  dplyr::full_join((mutationCalls %>%
                      dplyr::select(depmap_id, entrez_id,
                                    is_cosmic_hotspot, var_annotation)),
                   by = c("depmap_id", "entrez_id"))

p1 <- meta_crispr %>%
  ggplot(aes(x = dependency, y = lineage)) +
  geom_point(alpha = 0.4, size = 0.5) +
  geom_point(data = subset(
    meta_crispr, var_annotation == "damaging"), color = "red") +
  geom_point(data = subset(
    meta_crispr, var_annotation == "other non-conserving"), color = "blue") +
  geom_point(data = subset(
    meta_crispr, var_annotation == "other conserving"), color = "cyan") +
  geom_point(data = subset(
    meta_crispr, is_cosmic_hotspot == TRUE), color = "orange") +
  geom_vline(xintercept=mean(meta_crispr$dependency, na.rm = TRUE),
             linetype = "dotted", color = "red") +
  ggtitle("Scatterplot of dependency scores for gene CDKN2C by lineage")

p1

metadata %>%
  dplyr::select(depmap_id, lineage) %>%
  dplyr::full_join(TPM, by = "depmap_id") %>%
  dplyr::filter(gene_name == "MCM10") %>% 
  ggplot(aes(x = lineage, y = expression, fill = lineage)) +
  geom_boxplot(outlier.alpha = 0.1) +
  ggtitle("Boxplot of expression values for gene MCM10 by lineage") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.position = "none")

## expression vs rnai gene dependency for Glioblastoma
glioma <- metadata %>%
  dplyr::select(depmap_id, cell_line,
                primary_disease, subtype_disease) %>%
  dplyr::filter(primary_disease == "Brain Cancer",
                subtype_disease == "Glioblastoma")

rnai_sub <- rnai %>% dplyr::select(depmap_id, gene, gene_name, dependency)
tpm_sub <- TPM %>% dplyr::select(depmap_id, gene, gene_name, expression)
crispr_sub <- crispr %>% dplyr::select(depmap_id, gene, gene_name, dependency)

glioma_dep <- glioma %>%
  dplyr::left_join(crispr_sub, by = "depmap_id") %>%
  dplyr::select(-cell_line, -primary_disease,
                -subtype_disease, -gene_name)

glioma_exp <- glioma %>% dplyr::left_join(tpm_sub, by = "depmap_id")

glioma_dat_exp <- dplyr::full_join(glioma_dep, glioma_exp,
                                    by = c("depmap_id", "gene")) %>%
  dplyr::filter(!is.na(expression))

p2 <- ggplot(data = glioma_dat_exp, aes(x = dependency, y = expression)) +
  geom_point(alpha = 0.4, size = 0.5) +
  geom_vline(xintercept=mean(glioma_dat_exp$dependency, na.rm = TRUE),
             linetype = "dotted", color = "red") +
  geom_hline(yintercept=mean(glioma_dat_exp$expression, na.rm = TRUE),
             linetype = "dotted", color = "red") +
  ggtitle("Scatterplot of CRISPR dependency vs expression values for gene")
p2 + theme(axis.text.x = element_text(angle = 45))

glioma_dat_exp %>%
  dplyr::select(cell_line, gene_name, dependency, expression) %>%
  dplyr::arrange(desc(dependency), expression) %>% 
  head(10)

metadata %>%
  dplyr::select(depmap_id, lineage) %>%
  dplyr::full_join(copyNumber, by = "depmap_id") %>%
  dplyr::filter(gene_name == "MCM10") %>%
  ggplot(aes(x = lineage, y = log_copy_number, fill = lineage)) +
  geom_boxplot(outlier.alpha = 0.1) +
  ggtitle("Boxplot of log copy number for gene MCM10 by lineage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
