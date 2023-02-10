d <- read_csv("results/rna_results.csv") %>%
  select(entrez, log2FoldChange) %>%
  distinct(entrez, .keep_all = TRUE)
## assume 1st column is ID
## 2nd column is FC

## feature 1: numeric vector
geneList <- d$log2FoldChange

## feature 2: named vector
names(geneList) <- pull(d, entrez)

## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = TRUE)

goplot(ego)