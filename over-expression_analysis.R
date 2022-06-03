# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)

## Merge the AnnotationHub dataframe with the results 
res_ids <- left_join(res_tableHN_tb, annotations_ahb, by=c("gene"="gene_id_version"))

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allHN_genes <- as.character(res_ids$gene_id)

## Extract significant results
sigHN <- dplyr::filter(res_ids, padj < 0.05)

sigHN_genes <- as.character(sigHN$gene_id)

## Run GO enrichment analysis 
ego <- enrichGO(gene = sigHN_genes, 
                universe = allHN_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "results/clusterProfiler_HN.csv")

## Dotplot 
dotplot(ego, showCategory=50)

## Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50)


