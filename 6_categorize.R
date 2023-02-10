### STEP 6: Categorizing genes by transcription and translation regulation ###

# Load libraries
library(tidyverse)

merge_rnarpf <- dplyr::inner_join(
  x = resdf_rna,
  y = resdf_rpf,
  by = c("entrez", "symbol"),
  suffix = c(".rna", ".rpf")) %>%
  tidyr::drop_na(padj.rna, padj.rpf, log2FoldChange.rna, log2FoldChange.rpf)

merge_rnarpf <- merge_rnarpf %>% # strict p-value and LFC thresholds
  dplyr::mutate(category_rna = dplyr::case_when(
    padj.rna < 1e-05 & log2FoldChange.rna >= log2(2) ~ "H",
    padj.rna < 1e-05 & log2FoldChange.rna <= -log2(2) ~ "L",
    padj.rna < 1e-05 & abs(log2FoldChange.rna) < log2(2) ~ "M",
    padj.rna > 1e-05  ~ "M")
  ) %>%
  dplyr::mutate(category_rpf = dplyr::case_when(
    padj.rpf < 1e-05 & log2FoldChange.rpf >= log2(2) ~ "H",
    padj.rpf < 1e-05 & log2FoldChange.rpf <= -log2(2) ~ "L",
    padj.rpf < 1e-05 & abs(log2FoldChange.rpf) < log2(2) ~ "M",
    padj.rpf > 1e-05  ~ "M")
  ) %>%
  tidyr::unite(category_rnarpf, c("category_rna", "category_rpf"))

# Remove genes that do not fall into clusters
merge_rnarpf <- merge_rnarpf %>%
  dplyr::filter(!(category_rnarpf == "H_H" & (log2FoldChange.rna < log2(2) | log2FoldChange.rpf < log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "H_L" & (log2FoldChange.rna < log2(2) | log2FoldChange.rpf > -log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "H_M" & (log2FoldChange.rna < log2(2) | abs(log2FoldChange.rpf) > log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "L_H" & (log2FoldChange.rna > log2(2) | log2FoldChange.rpf < log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "L_L" & (log2FoldChange.rna > log2(2) | log2FoldChange.rpf > -log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "L_M" & (log2FoldChange.rna > log2(2) | abs(log2FoldChange.rpf) > log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "M_H" & (abs(log2FoldChange.rna) > log2(2) | log2FoldChange.rpf < log2(2)))) %>%
  dplyr::filter(!(category_rnarpf == "M_L" & (abs(log2FoldChange.rna) > log2(2) | log2FoldChange.rpf > -log2(2))))

merge_rnarpf <- left_join(merge_rnarpf, delta_te, by = c("entrez" = "gene_id")) %>%
  tidyr::drop_na(log2_delta_te, delta_te) %>%
  mutate(category_te = dplyr::case_when(
    category_rnarpf == "L_H" | category_rnarpf == "M_H"  ~ "TE.up",
    category_rnarpf == "H_L" | category_rnarpf == "M_L"  ~ "TE.down",
    category_rnarpf == "H_H" | category_rnarpf == "L_L"  ~ "TE.con",
    category_rnarpf == "L_M" | category_rnarpf == "H_M" ~ "TE.buff",
    category_rnarpf == "MM" ~ "NS"
    )
  )

# merge.te <- dplyr::inner_join(x = resdf.rna,
#                             y = resdf.rpf,
#                             by = "gene_id",
#                             suffix = c(".rna", ".rpf")) %>%
#   left_join(delta_te, by = "gene_id") %>%
#   drop_na(padj.rna,
#           padj.rpf,
#           log2FoldChange.rna,
#           log2FoldChange.rpf,
#           log2.delta_te,
#           delta_te) %>%
#   mutate(level.te = dplyr::case_when(
#     log2.delta_te >= log2(1) ~ "TE.up",
#     log2.delta_te <= log2(1) ~ "TE.down")
#     #log2.delta_te < log2(1.5) & log2.delta_te > log2(2/3) ~ "TE.mid")
#     )

write_csv(merge_rnarpf,
          file = "results/te_categories.csv")

scatterNums <- summarise(merge_rnarpf %>% group_by(category_rnarpf), num = n())
scatterNums

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


scatterlfc <- ggplot(merge_rnarpf %>% filter(category_rnarpf != "M_M"),
                     aes(x = log2FoldChange.rna,
                         y = log2FoldChange.rpf,
                         color = category_rnarpf)) +
  geom_point() +
  scale_colour_brewer(palette = "Dark2",
                      name = "Category",
                      labels = c("HH", "HL", "HM", "LH", "LL", "LM", "MH", "ML")) +
  labs(x = expression(Delta*"RNA"),
       y = expression(Delta*"RPF")) +
  theme_classic()
scatterlfc

# ggplot(scatterData %>% filter(category.rna.rpf != "M_M"),
#        aes(x=category.rna.rpf, y=log2.delta_te)) +
#   geom_boxplot(outlier.shape = NA) +
#   coord_cartesian(ylim=c(0,6)) +
#   labs(x = "Category",
#        y = expression(Delta*"TE")) +
#   scale_x_discrete(labels = c("HH", "HL", "HM", "LH", "LL", "LM", "MH", "ML")) +
#   geom_hline(aes(yintercept=1), linetype = "dashed", color = "red") +
#   theme_classic()

box <- ggplot(merge_rnarpf %>% filter(category_rnarpf != "M_M"),
       aes(x = category_rnarpf, y = log2_delta_te, fill = category_rnarpf)) +
  geom_boxplot() +
  coord_cartesian() +
  labs(x = "Category",
       y = expression(Delta*"TE")) +
  scale_fill_brewer(palette="Dark2") +
  scale_x_discrete(labels = c("HH", "HL", "HM", "LH", "LL", "LM", "MH", "ML")) +
  geom_hline(aes(yintercept=0), linetype = "dashed", color = "red") +
  theme_classic() +
  theme(legend.position = "none")
box

# tedata <- dplyr::select(te, gene_id, treatment, log2.te, te) %>%
#   pivot_wider(names_from = treatment, values_from = c(log2.te, te)) %>%
#   left_join(merge, by = "gene_id") %>%
#   drop_na(category.rna.rpf)
# 
# scatterte <- ggplot(tedata %>% filter(category.rna.rpf != "M_M"),
#                     aes(x=log2.te_Normoxic,
#                         y=log2.te_Hypoxic,
#                         color = category.rna.rpf)) +
#   geom_point() +
#   scale_colour_brewer(palette = "Dark2",
#                       name = "Category",
#                       labels = c("HH", "HL", "HM", "LH", "LL", "LM", "MH", "ML")) +
#   labs(x = "TE (Normoxia)",
#        y = "TE (Hypoxia)") +
#   theme_classic()
# scatterte

plot_grid(
  scatterlfc,
  pie,
  scatterte,
  box,
  axis = "l",
  align = 'vh',
  labels = "AUTO",
  nrow = 2
)

ggplot(merge_rnarpf %>% filter(!str_ends(.$category_rnarpf,"M")),
       aes(x=category_te, y=log2_delta_te)) +
  geom_boxplot() +
  coord_cartesian() +
  labs(x = "Category",
       y = expression(paste("log","2","(",Delta*"TE)"))) +
  #scale_x_discrete(labels = c("HH", "HL", "HM", "LH", "LL", "LM", "MH", "ML", "MM")) +
  geom_hline(aes(yintercept=0), linetype = "dashed", color = "red") +
  theme_classic()
