### STEP 7: Categorizing genes by transcription and translation regulation ###

## This script organizes genes into terciles based on levels of regulation in 
## RNA-seq and RIBO-seq

# Join RNA-seq and RIBO-seq results tables
merge_rnarpf_hn <- dplyr::inner_join(
  x = resdf_rna_hn,
  y = resdf_rpf_hn,
  by = c("entrez", "symbol"),
  suffix = c(".rna", ".rpf")) %>%
  tidyr::drop_na(padj.rna, padj.rpf, log2FoldChange.rna, log2FoldChange.rpf)

merge_rnarpf_gscnpc <- dplyr::inner_join(
  x = resdf_rna_gscnpc,
  y = resdf_rpf_gscnpc,
  by = c("entrez", "symbol"),
  suffix = c(".rna", ".rpf")) %>%
  tidyr::drop_na(padj.rna, padj.rpf, log2FoldChange.rna, log2FoldChange.rpf)

# Organize genes from RNA/RPF with one-letter codes and combine
merge_rnarpf_hn <- merge_rnarpf_hn %>% # strict p-value (0.10) and LFC (1) thresholds
  dplyr::mutate(category_rna = dplyr::case_when(
    padj.rna < 0.1 & log2FoldChange.rna >= log2(2) ~ "H",
    padj.rna < 0.1 & log2FoldChange.rna <= -log2(2) ~ "L",
    padj.rna < 0.1 & abs(log2FoldChange.rna) < log2(2) ~ "M",
    padj.rna > 0.1  ~ "M")
  ) %>%
  dplyr::mutate(category_rpf = dplyr::case_when(
    padj.rpf < 0.1 & log2FoldChange.rpf >= log2(2) ~ "H",
    padj.rpf < 0.1 & log2FoldChange.rpf <= -log2(2) ~ "L",
    padj.rpf < 0.1 & abs(log2FoldChange.rpf) < log2(2) ~ "M",
    padj.rpf > 0.1  ~ "M")
  ) %>%
  tidyr::unite(category_rnarpf, c("category_rna", "category_rpf"))

merge_rnarpf_gscnpc <- merge_rnarpf_gscnpc %>% # strict p-value (0.10) and LFC (1) thresholds
  dplyr::mutate(category_rna = dplyr::case_when(
    padj.rna < 0.1 & log2FoldChange.rna >= log2(2) ~ "H",
    padj.rna < 0.1 & log2FoldChange.rna <= -log2(2) ~ "L",
    padj.rna < 0.1 & abs(log2FoldChange.rna) < log2(2) ~ "M",
    padj.rna > 0.1  ~ "M")
  ) %>%
  dplyr::mutate(category_rpf = dplyr::case_when(
    padj.rpf < 0.1 & log2FoldChange.rpf >= log2(2) ~ "H",
    padj.rpf < 0.1 & log2FoldChange.rpf <= -log2(2) ~ "L",
    padj.rpf < 0.1 & abs(log2FoldChange.rpf) < log2(2) ~ "M",
    padj.rpf > 0.1  ~ "M")
  ) %>%
  tidyr::unite(category_rnarpf, c("category_rna", "category_rpf"))

# Remove genes that do not fall into clusters
## From here, only contintue with 'hypoxia versus normoxia' dataset
merge_rnarpf_hn <- merge_rnarpf_hn %>%
  dplyr::filter(!(category_rnarpf == "H_H" & (log2FoldChange.rna < log2(2) | log2FoldChange.rpf < log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "H_L" & (log2FoldChange.rna < log2(2) | log2FoldChange.rpf > -log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "H_M" & (log2FoldChange.rna < log2(2) | abs(log2FoldChange.rpf) > log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "L_H" & (log2FoldChange.rna > log2(2) | log2FoldChange.rpf < log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "L_L" & (log2FoldChange.rna > log2(2) | log2FoldChange.rpf > -log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "L_M" & (log2FoldChange.rna > log2(2) | abs(log2FoldChange.rpf) > log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "M_H" & (abs(log2FoldChange.rna) > log2(2) | log2FoldChange.rpf < log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "M_L" & (abs(log2FoldChange.rna) > log2(2) | log2FoldChange.rpf > -log2(2))))

merge_rnarpf_hn <- left_join(merge_rnarpf_hn, delta_te_hn, by = c("entrez" = "gene_id")) %>%
  tidyr::drop_na(log2_delta_te, delta_te) %>%
  mutate(category_te = dplyr::case_when(
    category_rnarpf == "L_H" | category_rnarpf == "M_H"  ~ "TE_up",
    category_rnarpf == "H_L" | category_rnarpf == "M_L"  ~ "TE_down",
    category_rnarpf == "H_H" | category_rnarpf == "L_L"  ~ "TE_con",
    category_rnarpf == "L_M" | category_rnarpf == "H_M" ~ "TE_buff",
    category_rnarpf == "M_M" ~ "NS"
    )
  )

write_csv(merge_rnarpf_hn,
          file = "results/te_hn_categories.csv")

scatterNums <- summarise(merge_rnarpf_hn %>% group_by(category_rnarpf), num = n())
scatterNums

# Visualize the proportions of genes in categories
pie <- ggplot(dplyr::slice(scatterNums, 1:8), aes(x = "", y = num, fill = category_rnarpf)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y", start = 0) +
  geom_label(aes(label = paste0(round(num/sum(dplyr::slice(scatterNums, 1:8)$num)*100,2), "%")),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
pie


scatterlfc <- ggplot(merge_rnarpf_hn %>% filter(category_rnarpf != "M_M"),
                     aes(x = log2FoldChange.rna,
                         y = log2FoldChange.rpf,
                         color = category_rnarpf)) +
  geom_point() +
  scale_colour_brewer(palette = "Dark2",
                      name = "Category",
                      labels = c("HH", "HL", "HM", "LH", "LL", "LM", "MH", "ML")) +
  labs(x = expression(Delta*"RNA (Log2 FC)"),
       y = expression(Delta*"RPF (Log2 FC)")) +
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -1, linetype = "dotted", color = "grey") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = -1, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey") +
  theme_classic()
scatterlfc

box <- ggplot(merge_rnarpf_hn %>% filter(category_rnarpf != "M_M"),
       aes(x = category_rnarpf, y = log2_delta_te, fill = category_rnarpf)) +
  geom_boxplot() +
  coord_cartesian() +
  labs(x = "Category",
       y = expression(Delta*"TE (Log2 FC Ratio)")) +
  scale_fill_brewer(palette="Dark2") +
  scale_x_discrete(labels = c("HH", "HL", "HM", "LH", "LL", "LM", "MH", "ML")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic() +
  theme(legend.position = "none")
box

plot_grid(
  scatterlfc,
  box,
  axis = "b",
  align = 'h',
  labels = "AUTO",
  nrow = 1
)
