### STEP 6: Running analysis and exploring DESeq2 results ###


# Perform differential expression testing with DESeq2 ---------------------
dds.rna_hn <- DESeq(dds.rna_hn)
dds.rpf_hn <- DESeq(dds.rpf_hn)

# Group cell lines into GSC or NPC categories for DESeq contrast
dds.rna_gscnpc$cell_type <- fct_collapse(dds.rna_gscnpc$cell_line,
                                         gsc = c("216", "315", "318"),
                                         npc = c("H1", "H9"))
dds.rna_gscnpc$cell_type <- relevel(dds.rna_gscnpc$cell_type, ref = "npc")
dds.rpf_gscnpc$cell_type <- fct_collapse(dds.rpf_gscnpc$cell_line,
                                         gsc = c("216", "315", "318"),
                                         npc = c("H1", "H9"))
dds.rpf_gscnpc$cell_type <- relevel(dds.rpf_gscnpc$cell_type, ref = "npc")
design(dds.rna_gscnpc) <- ~ cell_type
design(dds.rpf_gscnpc) <- ~ cell_type

dds.rna_gscnpc <- DESeq(dds.rna_gscnpc)
dds.rpf_gscnpc <- DESeq(dds.rpf_gscnpc)

# Build the results tables
rna_hn_res <- results(dds.rna_hn)
rpf_hn_res <- results(dds.rpf_hn)
rna_gscnpc_res <- results(dds.rna_gscnpc)
rpf_gscnpc_res <- results(dds.rpf_gscnpc)

mcols(rna_hn_res, use.names = TRUE) # information on each column in results
summary(rna_hn_res)
summary(rpf_hn_res)
summary(rna_gscnpc_res)
summary(rpf_gscnpc_res)

# Apply fold change shrinkage, replacing previous results table
resultsNames(dds.rna_hn)
rna_hn_res <- lfcShrink(dds = dds.rna_hn,
                        coef = "treatment_hypoxia_vs_normoxia",
                        res = rna_hn_res,
                        type = "apeglm")

resultsNames(dds.rpf_hn)
rpf_hn_res <- lfcShrink(dds = dds.rpf_hn,
                        coef = "treatment_hypoxia_vs_normoxia",
                        res = rpf_hn_res,
                        type = "apeglm")

resultsNames(dds.rna_gscnpc)
rna_gscnpc_res <- lfcShrink(dds = dds.rna_gscnpc,
                            coef = "cell_type_gsc_vs_npc",
                            res = rna_gscnpc_res,
                            type = "apeglm")

resultsNames(dds.rpf_gscnpc)
rpf_gscnpc_res <- lfcShrink(dds = dds.rpf_gscnpc,
                            coef = "cell_type_gsc_vs_npc",
                            res = rpf_gscnpc_res,
                            type = "apeglm")

# MA plot 
plotMA(rna_hn_res, ylim = c(-5, 5))
plotMA(rpf_hn_res, ylim = c(-5, 5))
plotMA(rna_gscnpc_res, ylim = c(-5, 5))
plotMA(rpf_gscnpc_res, ylim = c(-5, 5))

# Isolate the most significant gene in RIBO-seq 'hypoxia versus normoxia' contrast
# and display on the MA plot
rpf_hn_topGene <- rownames(rpf_hn_res)[which.min(rpf_hn_res$padj)]
with(rpf_hn_res[rpf_hn_topGene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, rpf_hn_topGene, pos = 2, col = "dodgerblue")
})

# Histogram of the p-values in RNA-seq 'hypoxia versus normoxia' contrast
hist(rna_hn_res$pvalue[rna_hn_res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

# Gene clustering using heatmap
topVarGenes_rna_hn <- head(order(rowVars(assay(vsd.rna_hn)), decreasing = TRUE), 20) # select top 20 genes rows with highest variance estimates
topVarGenes_rpf_hn <- head(order(rowVars(assay(vsd.rpf_hn)), decreasing = TRUE), 20)

mat_rna_hn  <- assay(vsd.rna_hn)[ topVarGenes_rna_hn, ]
mat_rna_hn  <- mat_rna_hn - rowMeans(mat_rna_hn)
anno_rna_hn <- as.data.frame(colData(vsd.rna_hn)[, c("cell_line","treatment")])
pheatmap(mat_rna_hn, annotation_col = anno_rna_hn) # heatmap of relative VST-transformed values across samples

mat_rpf_hn  <- assay(vsd.rpf_hn)[ topVarGenes_rpf_hn, ]
mat_rpf_hn  <- mat_rpf_hn - rowMeans(mat_rpf_hn)
anno_rpf_hn <- as.data.frame(colData(vsd.rpf_hn)[, c("cell_line","treatment")])
pheatmap(mat_rpf_hn, annotation_col = anno_rpf_hn)

# Convert normalized_counts to a data frame and transfer the row names to a column
normCounts_rna_hn <- counts(dds.rna_hn, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_id")

normCounts_rpf_hn <- counts(dds.rpf_hn, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_id")

normCounts_rna_gscnpc <- counts(dds.rna_gscnpc, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_id")

normCounts_rpf_gscnpc <- counts(dds.rpf_gscnpc, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_id")

save(dds.rna_hn, dds.rpf_hn,
     dds.rna_gscnpc, dds.rpf_gscnpc,
     rna_hn_res, rpf_hn_res,
     rna_gscnpc_res, rpf_gscnpc_res,
     normCounts_rna_hn, normCounts_rpf_hn,
     normCounts_rna_gscnpc, normCounts_rpf_gscnpc,
     vsd.rna_hn, vsd.rpf_hn,
     vsd.rna_gscnpc, vsd.rpf_gscnpc,
     all_hn_tpm, all_gscnpc_tpm,
     delta_te_hn,
     file = "data/RNA-RIBO_DESeqResults.RData")
