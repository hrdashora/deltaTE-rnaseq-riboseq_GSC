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
dotplot(ego, showCategory=50, font.size = 10)

## Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 30)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
HN_foldchanges <- sigHN$log2FoldChange

names(HN_foldchanges) <- sigHN$gene_id

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=HN_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
HN_foldchanges <- ifelse(HN_foldchanges > 2, 2, HN_foldchanges)
HN_foldchanges <- ifelse(HN_foldchanges < -2, -2, HN_foldchanges)

cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=HN_foldchanges, 
         vertex.label.font=6)



