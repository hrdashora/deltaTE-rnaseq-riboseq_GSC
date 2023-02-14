### STEP 4: Annotating significant genes ###
load("~/opt/r/rnaseq-riboseq/RNA-RIBO_SummarizedExperiment.RData")
load("~/opt/r/rnaseq-riboseq/RNA-RIBO_DESeqResults.RData")

# Load libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)

# Provide rownames of results table as a key
ens_str_rna <- rownames(rna_res)
ens_str_rpf <- rownames(rpf_res)

# Add individual columns to results table, ensure 'select()' performed 1-to-1 mapping
rna_res$symbol <- mapIds(org.Hs.eg.db,
                         keys = ens_str_rna,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
rpf_res$symbol <- mapIds(org.Hs.eg.db,
                         keys = ens_str_rpf,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")

resOrdered_rna <- rna_res[order(rna_res$pvalue),]
resOrdered_rpf <- rpf_res[order(rpf_res$pvalue),]

# Determine how many of the symbol column entries are NA
length(which(is.na(resOrdered_rna$symbol))) # 117 
length(which(is.na(resOrdered_rpf$symbol))) # 127

keep_rna <- !is.na(resOrdered_rna$symbol)
keep_rpf <- !is.na(resOrdered_rpf$symbol)

resKeepOrdered_rna <- resOrdered_rna[keep_rna, ]
resKeepOrdered_rpf <- resOrdered_rpf[keep_rpf, ]

# Create a tibble of results
resdf_rna <- resKeepOrdered_rna %>%
  as.data.frame() %>%
  rownames_to_column(var = "entrez") %>%
  as_tibble()
resdf_rpf <- resKeepOrdered_rpf %>%
  as.data.frame() %>%
  rownames_to_column(var = "entrez") %>%
  as_tibble()

write_excel_csv(resdf_rna, file = "results/rna_results.csv")
write_excel_csv(resdf_rpf, file = "results/rpf_results.csv")

save(rna_dds, rpf_dds,
     rna_res, rpf_res,
     resKeepOrdered_rna, resKeepOrdered_rpf,
     resdf_rna, resdf_rpf,
     normCounts_rna, normCounts_rpf,
     delta_te,
     file = "RNA-RIBO_DESeqResults.RData")

