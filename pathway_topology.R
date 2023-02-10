### STEP 8: Pathway Topology ###

library(SPIA)

# Significant genes is a vector of fold changes where the names are ENTREZ gene IDs.
# The background set is a vector of all the genes represented on the platform.

background_entrez.rna <- res_entrez.rna$gene_id
sig_res_entrez.rna <- res_entrez.rna[which(res_entrez.rna$padj < 0.05), ]
sig_entrez.rna <- sig_res_entrez.rna$log2FoldChange
names(sig_entrez.rna) <- sig_res_entrez.rna$gene_id

head(sig_entrez.rna)

spia_result <- spia(de=sig_entrez.rna, all=background_entrez.rna, organism="hsa")

head(spia_result, n=20)

plotP(spia_result, threshold=0.05)

## Look at specific pathway and view kegglink
subset(spia_result, ID == "04510")


#####


