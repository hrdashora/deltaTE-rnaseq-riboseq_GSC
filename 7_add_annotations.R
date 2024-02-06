### STEP 5: Annotating Entrez gene IDs with gene symbols ###

## This script will add gene symbols to the results tables and save them as CSV
## files and update objects in .RData files.

load("data/RNA-RIBO_SummarizedExperiment.RData")
load("data/RNA-RIBO_DESeqResults.RData")

# Load libraries
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(tidyverse)
})

# Provide rownames of results table as a key
ens_str_rna_hn <- rownames(rna_hn_res)
ens_str_rpf_hn <- rownames(rpf_hn_res)
ens_str_rna_gscnpc <- rownames(rna_gscnpc_res)
ens_str_rpf_gscnpc <- rownames(rpf_gscnpc_res)

# Add individual columns to results table, ensure 'select()' performed
# 1-to-1 mapping
rna_hn_res$symbol <- mapIds(org.Hs.eg.db,
                            keys = ens_str_rna_hn,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
rpf_hn_res$symbol <- mapIds(org.Hs.eg.db,
                            keys = ens_str_rpf_hn,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
rna_gscnpc_res$symbol <- mapIds(org.Hs.eg.db,
                                keys = ens_str_rna_gscnpc,
                                column = "SYMBOL",
                                keytype = "ENTREZID",
                                multiVals = "first")
rpf_gscnpc_res$symbol <- mapIds(org.Hs.eg.db,
                                keys = ens_str_rpf_gscnpc,
                                column = "SYMBOL",
                                keytype = "ENTREZID",
                                multiVals = "first")

resOrdered_rna_hn <- rna_hn_res[order(rna_hn_res$padj),] # order results by adjusted p-value
resOrdered_rpf_hn <- rpf_hn_res[order(rpf_hn_res$padj),]
resOrdered_rna_gscnpc <- rna_gscnpc_res[order(rna_gscnpc_res$padj),]
resOrdered_rpf_gscnpc <- rpf_gscnpc_res[order(rpf_gscnpc_res$padj),]

# Determine how many of results table values are NA
## Row contains a sample with an extreme count outlier
length(which(is.na(resOrdered_rna_hn$pvalue))) # 0
length(which(is.na(resOrdered_rpf_hn$pvalue))) # 0
length(which(is.na(resOrdered_rna_gscnpc$pvalue))) # 0
length(which(is.na(resOrdered_rpf_gscnpc$pvalue))) # 66

## Row is marked by automatic independent filtering for having a low mean normalized count
length(which(is.na(resOrdered_rna_hn$padj))) # 0
length(which(is.na(resOrdered_rpf_hn$padj))) # 0
length(which(is.na(resOrdered_rna_gscnpc$padj))) # 0
length(which(is.na(resOrdered_rpf_gscnpc$padj))) # 66

# Keep genes witch pass the independent filtering threshold
keep_rna_hn <- !is.na(resOrdered_rna_hn$padj)
keep_rpf_hn <- !is.na(resOrdered_rpf_hn$padj)
resKeepOrdered_rna_hn <- resOrdered_rna_hn[keep_rna_hn, ]
resKeepOrdered_rpf_hn <- resOrdered_rpf_hn[keep_rpf_hn, ]

keep_rna_gscnpc <- !is.na(resOrdered_rna_gscnpc$padj)
keep_rpf_gscnpc <- !is.na(resOrdered_rpf_gscnpc$padj)
resKeepOrdered_rna_gscnpc <- resOrdered_rna_gscnpc[keep_rna_gscnpc, ]
resKeepOrdered_rpf_gscnpc <- resOrdered_rpf_gscnpc[keep_rpf_gscnpc, ]

# How many of the symbol column entries are NA?
length(which(is.na(resKeepOrdered_rna_hn$symbol))) # 50
length(which(is.na(resKeepOrdered_rpf_hn$symbol))) # 75
length(which(is.na(resKeepOrdered_rna_gscnpc$symbol))) # 70
length(which(is.na(resKeepOrdered_rpf_gscnpc$symbol))) # 114

# Create a tibble of results
resdf_rna_hn <- resKeepOrdered_rna_hn %>%
  as.data.frame() %>%
  rownames_to_column(var = "entrez") %>%
  as_tibble()
resdf_rpf_hn <- resKeepOrdered_rpf_hn %>%
  as.data.frame() %>%
  rownames_to_column(var = "entrez") %>%
  as_tibble()
resdf_rna_gscnpc <- resKeepOrdered_rna_gscnpc %>%
  as.data.frame() %>%
  rownames_to_column(var = "entrez") %>%
  as_tibble()
resdf_rpf_gscnpc <- resKeepOrdered_rpf_gscnpc %>%
  as.data.frame() %>%
  rownames_to_column(var = "entrez") %>%
  as_tibble()

write_excel_csv(resdf_rna_hn, file = "results/rna_hn_results.csv")
write_excel_csv(resdf_rpf_hn, file = "results/rpf_hn_results.csv")
write_excel_csv(resdf_rna_gscnpc, file = "results/rna_gscnpc_results.csv")
write_excel_csv(resdf_rpf_gscnpc, file = "results/rpf_gscnpc_results.csv")

# Re-save important objects into .RData file
save(dds.rna_hn, dds.rpf_hn,
     dds.rna_gscnpc, dds.rpf_gscnpc,
     rna_hn_res, rpf_hn_res,
     rna_gscnpc_res, rpf_gscnpc_res,
     resKeepOrdered_rna_hn, resKeepOrdered_rpf_hn,
     resKeepOrdered_rna_gscnpc, resKeepOrdered_rpf_gscnpc,
     normCounts_rna_hn, normCounts_rpf_hn,
     normCounts_rna_gscnpc, normCounts_rpf_gscnpc,
     vsd.rna_hn, vsd.rpf_hn,
     vsd.rna_gscnpc, vsd.rpf_gscnpc,
     all_hn_tpm, all_gscnpc_tpm,
     delta_te_hn,
     file = "data/RNA-RIBO_DESeqResults.RData")
