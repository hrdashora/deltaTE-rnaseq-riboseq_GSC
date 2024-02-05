### STEP 5: Annotating Entrez gene IDs with gene symbols ###

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

resOrdered_rna <- rna_res[order(rna_res$padj),] # order results by adjusted p-value
resOrdered_rpf <- rpf_res[order(rpf_res$padj),]

# Determine how many of results table values are NA
## Row contains a sample with an extreme count outlier
length(which(is.na(resOrdered_rna$pvalue))) # 0
length(which(is.na(resOrdered_rpf$pvalue)))

## Row is marked by automatic independent filtering for having a low mean normalized count
length(which(is.na(resOrdered_rna$padj))) # 4513
length(which(is.na(resOrdered_rpf$padj))) # 2190

# Keep genes witch pass the independent filtering threshold
keep_rna <- !is.na(resOrdered_rna$padj)
keep_rpf <- !is.na(resOrdered_rpf$padj)
resKeepOrdered_rna <- resOrdered_rna[keep_rna, ]
resKeepOrdered_rpf <- resOrdered_rpf[keep_rpf, ]

# How many of the symbol column entries are NA?
length(which(is.na(resKeepOrdered_rna$symbol))) # 85
length(which(is.na(resKeepOrdered_rpf$symbol))) # 106

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

