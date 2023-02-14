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
