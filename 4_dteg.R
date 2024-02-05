### STEP 4: Calculate deltaTE using the DTEG.R framework ###

## This script will utilize the protocol outlined in DTEG.R to calculate deltaTE
## and complete the appropriate visualizations. All cell lines should be
## analyzed separately.


# Design DESeq to calculate DTEGs -----------------------------------------
# Prepare input files by subsetting SummarizedExperiment objects
rna_counts_hn <- assay(subset(se.rnaseq,
                              select = source == "JSY"))
colnames(rna_counts_hn) <- colData(subset(se.rnaseq,
                                          select = source == "JSY"))$file_name

ribo_counts_hn <- assay(subset(se.riboseq,
                               select = source == "JSY"))
colnames(ribo_counts_hn) <- colData(subset(se.riboseq,
                                           select = source == "JSY"))$file_name

rna_counts_gscnpc <- assay(subset(se.rnaseq,
                                  select = treatment == "normoxia")[,3:15]) # remove 318 replicate 1 and 2 from further analysis
colnames(rna_counts_gscnpc) <- colData(subset(se.rnaseq,
                                              select = treatment == "normoxia")[,3:15])$file_name

ribo_counts_gscnpc <- assay(subset(se.riboseq,
                                   select = treatment == "normoxia"))
colnames(ribo_counts_gscnpc) <- colData(subset(se.riboseq,
                                               select = treatment == "normoxia"))$file_name

# Create sample information data frames
colData(se.rnaseq)$assay <- factor(rep("rna", ncol(se.rnaseq)))
colData(se.riboseq)$assay <- factor(rep("ribo", ncol(se.riboseq)))

sample_info_hn <- data.frame(
  rbind(colData(subset(se.rnaseq,
                       select = source == "JSY")),
        colData(subset(se.riboseq,
                       select = source == "JSY")))) %>%
  dplyr::relocate(SampleID = file_name, Condition = treatment, SeqType = assay, Batch = cell_line) %>%
  dplyr::select(-c(replicate, source))

sample_info_gscnpc <- data.frame(
  rbind(colData(subset(se.rnaseq,
                       select = treatment == "normoxia")[,3:15]),
        colData(subset(se.riboseq,
                       select = treatment == "normoxia")))) %>%
  mutate(cell_type = fct_collapse(cell_line,
                                  gsc = c("216", "315", "318"),
                                  npc = c("H1", "H9"))) %>% # collapse cell lines into cell types
  dplyr::relocate(SampleID = file_name, Condition = cell_type, SeqType = assay, Batch = cell_line) %>%
  dplyr::select(-c(replicate, treatment, source))
sample_info_gscnpc$Condition <- fct_relevel(sample_info_gscnpc$Condition, "npc", "gsc")

# Merge counts data and create DESeq2 objects
identical(rownames(rna_counts_hn), rownames(ribo_counts_hn))
identical(sample_info_hn$SampleID, colnames(cbind(rna_counts_hn, ribo_counts_hn)))
ddsMat_hn <- DESeqDataSetFromMatrix(
  countData = cbind(rna_counts_hn, ribo_counts_hn),
  colData = sample_info_hn,
  design = ~ Batch + Condition + SeqType + Condition:SeqType
)

identical(rownames(rna_counts_gscnpc), rownames(ribo_counts_gscnpc))
identical(sample_info_gscnpc$SampleID, colnames(cbind(rna_counts_gscnpc, ribo_counts_gscnpc)))
ddsMat_gscnpc <- DESeqDataSetFromMatrix(
  countData = cbind(rna_counts_gscnpc, ribo_counts_gscnpc),
  colData = sample_info_gscnpc,
  design = ~ Condition + SeqType + Condition:SeqType
)

# Run DESeq2
ddsMat_hn <- DESeq(ddsMat_hn)
ddsMat_gscnpc <- DESeq(ddsMat_gscnpc)

# Obtain fold changes for TE
## Interaction term quantifies the change in TE in treatment versus control
resultsNames(ddsMat_hn)
res_hn <- results(ddsMat_hn, name = "Conditionhypoxia.SeqTyperibo")

resultsNames(ddsMat_gscnpc)
res_gscnpc <- results(ddsMat_gscnpc, name = "Conditiongsc.SeqTyperibo")

# Store the list of DTEGs in a file
write.table(res_hn[which(res_hn$padj < 0.05), ],
            "results/DTEGs_hn.txt",
            quote = FALSE)

write.table(res_gscnpc[which(res_gscnpc$padj < 0.05), ],
            "results/DTEGs_gscnpc.txt",
            quote = FALSE)

# Run DESeq2 for mRNA counts in order to obtain DTGs
ddsMat_hn_rna <- DESeqDataSetFromMatrix(
  countData = rna_counts_hn,
  colData = sample_info_hn[which(sample_info_hn$SeqType == "rna"),],
  design = ~ Batch + Condition
)
ddsMat_hn_rna <- DESeq(ddsMat_hn_rna)
resultsNames(ddsMat_hn_rna)
res_hn_rna <- lfcShrink(ddsMat_hn_rna,
                        coef = "Condition_hypoxia_vs_normoxia",
                        type = "apeglm")
write.table(res_hn_rna[which(res_hn_rna$padj < 0.05), ],
            "results/DEGs_hn.txt",
            quote = FALSE)

ddsMat_gscnpc_rna <- DESeqDataSetFromMatrix(
  countData = rna_counts_gscnpc,
  colData = sample_info_gscnpc[which(sample_info_gscnpc$SeqType == "rna"),],
  design = ~ Condition
)
ddsMat_gscnpc_rna <- DESeq(ddsMat_gscnpc_rna)
resultsNames(ddsMat_gscnpc_rna)
res_gscnpc_rna <- lfcShrink(ddsMat_gscnpc_rna,
                            coef = "Condition_gsc_vs_npc",
                            type = "apeglm")
write.table(res_gscnpc_rna[which(res_gscnpc_rna$padj < 0.05), ],
            "results/DEGs_gscnpc.txt",
            quote = FALSE)

# Run DESeq2 for RPFs
ddsMat_hn_ribo <- DESeqDataSetFromMatrix(
  countData = ribo_counts_hn,
  colData = sample_info_hn[which(sample_info_hn$SeqType == "ribo"),],
  design = ~ Batch + Condition
)
ddsMat_hn_ribo <- DESeq(ddsMat_hn_ribo)
resultsNames(ddsMat_hn_ribo)
res_hn_ribo <- lfcShrink(ddsMat_hn_ribo,
                         coef = "Condition_hypoxia_vs_normoxia",
                         type = "apeglm")
write.table(res_hn_ribo[which(res_hn_ribo$padj < 0.05), ],
            "results/RPFs_hn_ribo.txt",
            quote = FALSE)

ddsMat_gscnpc_ribo <- DESeqDataSetFromMatrix(
  countData = ribo_counts_gscnpc,
  colData = sample_info_gscnpc[which(sample_info_gscnpc$SeqType == "ribo"),],
  design = ~ Condition
)
ddsMat_gscnpc_ribo <- DESeq(ddsMat_gscnpc_ribo)
resultsNames(ddsMat_gscnpc_ribo)
res_gscnpc_ribo <- lfcShrink(ddsMat_gscnpc_ribo,
                             coef = "Condition_gsc_vs_npc",
                             type = "apeglm")
write.table(res_gscnpc_ribo[which(res_gscnpc_ribo$padj < 0.05), ],
            "results/RPFs_gscnpc_ribo.txt",
            quote = FALSE)

# Categorize genes into distinct regulation groups ------------------------
## Forwarded: genes driven by transcriptional regulation, change in RPF is in the same direction as change in RNA
forwarded_gscnpc <- rownames(res_gscnpc)[which(res_gscnpc$padj > 0.05 &
                                                 res_gscnpc_ribo$padj < 0.05 &
                                                 res_gscnpc_rna$padj < 0.05)]
forwarded_hn <- rownames(res_hn)[which(res_hn$padj > 0.05 &
                                         res_hn_ribo$padj < 0.05 &
                                         res_hn_rna$padj < 0.05)]
## Exclusive: genes regulated exclusively by translation, change in RPF is not driven by change in RNA
exclusive_gscnpc <- rownames(res_gscnpc)[which(res_gscnpc$padj < 0.05 & 
                                                 res_gscnpc_ribo$padj < 0.05 & 
                                                 res_gscnpc_rna$padj > 0.05)]
exclusive_gscnpc_pos <- rownames(res_gscnpc)[which(res_gscnpc$log2FoldChange > 0 &
                                                     res_gscnpc$padj < 0.05 & 
                                                     res_gscnpc_ribo$padj < 0.05 & 
                                                     res_gscnpc_rna$padj > 0.05)]
exclusive_gscnpc_neg <- rownames(res_gscnpc)[which(res_gscnpc$log2FoldChange < 0 &
                                                     res_gscnpc$padj < 0.05 & 
                                                     res_gscnpc_ribo$padj < 0.05 & 
                                                     res_gscnpc_rna$padj > 0.05)]

exclusive_hn <- rownames(res_hn)[which(res_hn$padj < 0.05 &
                                         res_hn_ribo$padj < 0.05 &
                                         res_hn_rna$padj > 0.05)]
exclusive_hn_pos <- rownames(res_hn)[which(res_hn$log2FoldChange > 0 &
                                             res_hn$padj < 0.05 &
                                             res_hn_ribo$padj < 0.05 &
                                             res_hn_rna$padj > 0.05)]
exclusive_hn_neg <- rownames(res_hn)[which(res_hn$log2FoldChange > 0 &
                                             res_hn$padj < 0.05 &
                                             res_hn_ribo$padj < 0.05 &
                                             res_hn_rna$padj > 0.05)]
## Intensified & Buffered: genes regulated by both transcriptional and translational regulation
both_gscnpc <- rownames(res_gscnpc)[which(res_gscnpc$padj < 0.05 &
                                            res_gscnpc_ribo$padj < 0.05 &
                                            res_gscnpc_rna$padj < 0.05)]
both_gscnpc_pos <- rownames(res_gscnpc)[which(res_gscnpc$log2FoldChange > 0 &
                                                res_gscnpc$padj < 0.05 &
                                                res_gscnpc_ribo$padj < 0.05 &
                                                res_gscnpc_rna$padj < 0.05)]
both_gscnpc_neg <- rownames(res_gscnpc)[which(res_gscnpc$log2FoldChange < 0 &
                                                res_gscnpc$padj < 0.05 &
                                                res_gscnpc_ribo$padj < 0.05 &
                                                res_gscnpc_rna$padj < 0.05)]
intensified_gscnpc <- both_gscnpc[which(res_gscnpc[both_gscnpc,2]*res_gscnpc_rna[both_gscnpc,2] > 0)]
intensified_gscnpc_pos <- both_gscnpc_pos[which(res_gscnpc[both_gscnpc_pos,2]*res_gscnpc_rna[both_gscnpc_pos,2] > 0)]
intensified_gscnpc_neg <- both_gscnpc_neg[which(res_gscnpc[both_gscnpc_neg,2]*res_gscnpc_rna[both_gscnpc_neg,2] > 0)]

buffered_gscnpc <- both_gscnpc[which(res_gscnpc[both_gscnpc,2]*res_gscnpc_rna[both_gscnpc,2] < 0)]
buffered_gscnpc <- c(rownames(res_gscnpc)[which(res_gscnpc$padj < 0.05 &
                                                  res_gscnpc_ribo$padj > 0.05 &
                                                  res_gscnpc_rna$padj < 0.05)],
                     buffered_gscnpc)

buffered_gscnpc_pos <- both_gscnpc_pos[which(res_gscnpc[both_gscnpc_pos,2]*res_gscnpc_rna[both_gscnpc_pos,2] < 0)]
buffered_gscnpc_pos <- c(rownames(res_gscnpc)[which(res_gscnpc$log2FoldChange > 0 &
                                                      res_gscnpc$padj < 0.05 &
                                                      res_gscnpc_ribo$padj > 0.05 &
                                                      res_gscnpc_rna$padj < 0.05)],
                         buffered_gscnpc_pos)

buffered_gscnpc_neg <- both_gscnpc_neg[which(res_gscnpc[both_gscnpc_neg,2]*res_gscnpc_rna[both_gscnpc_neg,2] < 0)]
buffered_gscnpc_neg <- c(rownames(res_gscnpc)[which(res_gscnpc$log2FoldChange < 0 &
                                                      res_gscnpc$padj < 0.05 &
                                                      res_gscnpc_ribo$padj > 0.05 &
                                                      res_gscnpc_rna$padj < 0.05)],
                         buffered_gscnpc_neg)

both_hn <- rownames(res_hn)[which(res_hn$padj < 0.05 &
                                    res_hn_ribo$padj < 0.05 &
                                    res_hn_rna$padj < 0.05)]
intensified_hn <- both_hn[which(res_hn[both_hn,2]*res_hn_rna[both_hn,2] > 0)]
buffered_hn <- both_hn[which(res_hn[both_hn,2]*res_hn_rna[both_hn,2] < 0)]
buffered_hn <- c(rownames(res_hn)[which(res_hn$padj < 0.05 &
                                          res_hn_ribo$padj > 0.05 &
                                          res_hn_rna$padj < 0.05)],
                 buffered_hn)

# Visualize global translation and transcriptional regulation
max_val <- max(res_gscnpc_ribo[,2], res_gscnpc_rna[,2], na.rm = T)
plot(y = res_gscnpc_ribo[,2],
     x = res_gscnpc_rna[,2],
     xlab = "RNA-seq log2 fold change",
     ylab = "Ribo-seq log2 fold change",
     asp = 1, 
     pch = 16,
     col = rgb(128/255, 128/255, 128/255, 0.1),
     ylim = c(-max_val, max_val),
     xlim = c(-max_val, max_val),
     cex = 0.4
)
abline(a = 0, b = 1, col = "gray")
abline(h = 0, v = 0, col = "gray")
points(y = res_gscnpc_ribo[forwarded_gscnpc, 2],
       x = res_gscnpc_rna[forwarded_gscnpc, 2],
       pch = 16,
       col = rgb(0,0,1,1),
       cex = 0.4)
points(y = res_gscnpc_ribo[exclusive_gscnpc, 2],
       x = res_gscnpc_rna[exclusive_gscnpc, 2],
       pch = 16, 
       col = rgb(1,0,0,1),
       cex = 0.4)
points(y = res_gscnpc_ribo[intensified_gscnpc, 2],
       x = res_gscnpc_rna[intensified_gscnpc, 2],
       pch = 16, 
       col = rgb(1,1,0,1),
       cex = 0.4)
points(y = res_gscnpc_ribo[buffered_gscnpc, 2],
       x = res_gscnpc_rna[buffered_gscnpc, 2],
       pch = 16, 
       col = rgb(1,0,1,1),
       cex = 0.4)

max_val <- max(res_hn_ribo[,2], res_hn_rna[,2], na.rm = T)
plot(y = res_hn_ribo[,2],
     x = res_hn_rna[,2],
     xlab = "RNA-seq log2 fold change",
     ylab = "Ribo-seq log2 fold change",
     asp = 1, 
     pch = 16,
     col = rgb(128/255, 128/255, 128/255, 0.1),
     ylim = c(-max_val, max_val),
     xlim = c(-max_val, max_val),
     cex = 0.4
)
abline(a = 0, b = 1, col = "gray")
abline(h = 0, v = 0, col = "gray")
points(y = res_hn_ribo[forwarded_hn, 2],
       x = res_hn_rna[forwarded_hn, 2],
       pch = 16, 
       col = rgb(0,0,1,1),
       cex = 0.4)
points(y = res_hn_ribo[exclusive_hn, 2],
       x = res_hn_rna[exclusive_hn, 2],
       pch = 16, 
       col = rgb(1,0,0,1),
       cex = 0.4)
points(y = res_hn_ribo[intensified_hn, 2],
       x = res_hn_rna[intensified_hn, 2],
       pch = 16, 
       col = rgb(1,1,0,1),
       cex = 0.4)
points(y = res_hn_ribo[buffered_hn, 2],
       x = res_hn_rna[buffered_hn, 2],
       pch = 16, 
       col = rgb(1,0,1,1),
       cex = 0.4)

