### STEP 3: Exploring DESeq2 results ###

# Perform differential expression testing with DESeq2
rna_dds <- DESeq(dds.rna)
rpf_dds <- DESeq(dds.rpf)

# Build the results table
rna_res <- results(rna_dds)
rpf_res <- results(rpf_dds)

mcols(rna_res, use.names = TRUE) # information on each column in results
summary(rna_res)
summary(rpf_res)

rna_res.05 <- results(rna_dds, alpha = 0.05)
rpf_res.05 <- results(rpf_dds, alpha = 0.05)
table(rna_res.05$padj < 0.05)

rna_resLFC1 <- results(rna_dds, lfcThreshold = 1)
rpf_resLFC1 <- results(rpf_dds, lfcThreshold = 1)
table(rna_resLFC1$padj < 0.1)

# Multiple testing
rna_resSig <- subset(rna_res, padj < 0.1)
head(rna_resSig[ order(rna_resSig$log2FoldChange), ]) # significant genes with the strongest down-regulation
head(rna_resSig[ order(rna_resSig$log2FoldChange, decreasing = TRUE), ]) # significant genes with the strongest up-regulation
rpf_resSig <- subset(rpf_res, padj < 0.1)
head(rpf_resSig[ order(rpf_resSig$log2FoldChange), ])
head(rpf_resSig[ order(rpf_resSig$log2FoldChange, decreasing = TRUE), ])

# Apply fold change shrinkage
rna_res <- lfcShrink(rna_dds,
                     coef = "treatment_Hypoxic_vs_Normoxic",
                     type = "apeglm")
rpf_res <- lfcShrink(rpf_dds,
                     coef = "treatment_Hypoxic_vs_Normoxic",
                     type = "apeglm")

# MA plot 
plotMA(rna_res, ylim = c(-5, 5))
plotMA(rpf_res, ylim = c(-5, 5))
topGene <- rownames(rpf_res)[which.min(rpf_res$padj)]
with(rpf_res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# Histogram of the p-values
hist(rna_res$pvalue[rna_res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

# Gene clustering
topVarGenes_rna <- head(order(rowVars(assay(vsd.rna)), decreasing = TRUE), 20)
topVarGenes_rpf <- head(order(rowVars(assay(vsd.rpf)), decreasing = TRUE), 20)

mat_rna  <- assay(vsd.rna)[ topVarGenes_rna, ]
mat_rna  <- mat_rna - rowMeans(mat_rna)
anno_rna <- as.data.frame(colData(vsd.rna)[, c("cell_line","treatment")])
pheatmap(mat_rna, annotation_col = anno_rna) # heatmap of relative VST-transformed values across samples

mat_rpf  <- assay(vsd.rpf)[ topVarGenes_rpf, ]
mat_rpf  <- mat_rpf - rowMeans(mat_rpf)
anno_rpf <- as.data.frame(colData(vsd.rpf)[, c("cell_line","treatment")])
pheatmap(mat_rpf, annotation_col = anno_rpf)

# Convert normalized_counts to a data frame and transfer the row names to a column
normCounts_rna <- counts(rna_dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_id")

normCounts_rpf <- counts(rpf_dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene_id")

save(rna_dds, rpf_dds,
     rna_res, rpf_res,
     normCounts_rna, normCounts_rpf,
     delta_te, file = "RNA-RIBO_DESeqResults.RData")
