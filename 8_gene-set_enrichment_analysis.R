### STEP 8: Gene Set Enrichment Analysis ###

# Create ordered ranked gene lists
sum(is.na(resTested$gene_id)) # Check if NA values are removed
sum(duplicated(resTested$gene_id) == T) # Check for Entrez ID duplicates

sum(is.na(sig_genes.rna$entrez))
sum(duplicated(is.na(sig_genes.rna$entrez) == T))

# Extract the fold changes and sort in decreasing order
foldchanges.rna <- resTested %>%
  dplyr::select(entrez, log2FoldChange.rna) %>%
  arrange(desc(log2FoldChange.rna)) %>%
  deframe()
foldchanges.rpf <- resTested %>%
  dplyr::select(entrez, log2FoldChange.rpf) %>%
  arrange(desc(log2FoldChange.rpf)) %>%
  deframe()

categoryRank <- sapply(categoryList, simplify = FALSE, USE.NAMES = TRUE, function(df){
  df %>%
    dplyr::select(gene_id, log2FoldChange.rpf) %>% # gene lists will be ordered by RPF in the category analysis
    arrange(desc(log2FoldChange.rpf)) %>%
    deframe() # convert two-column data frame into a named vector
})

set.seed(123456)

# KEGG pathway gene set enrichment analysis ---------------------------------------------------------------

# kegglist <- vector(mode='list', length(catgenelists))
# names(kegglist) <- names(catgenelists)
# keggsummary <- vector(mode='list', length(catgenelists))
# names(keggsummary) <- names(catgenelists)
# 
# # Create for loop to conduct GO analysis on categories and output results 
# # to a data frame
# for (j in 1:length(catgenelists)) {
#   input <- catgenelists[[j]]
#   kegglist[[j]] <- gseKEGG(gene = as.character(input$gene_id),
#                            universe = all_genes,
#                            keyType = "ENTREZID",
#                            OrgDb = org.Hs.eg.db,
#                            ont = "BP", # biological process
#                            pAdjustMethod = "BH",
#                            qvalueCutoff = 0.05,
#                            readable = TRUE)
#   
#   keggsummary[[j]] <- data.frame(kegglist[[j]])
# }

# GSEA using gene sets from KEGG pathways
gseaKEGG.rna <- gseKEGG(geneList = foldchanges.rna, # ordered named vector of fold changes (Entrez IDs are the associated names)
                        organism = "hsa", # supported organisms listed below
                        minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                        pvalueCutoff = 0.05, # padj cutoff value
                        verbose = TRUE)
gseaKEGG.rpf <- gseKEGG(geneList = foldchanges.rpf, # ordered named vector of fold changes (Entrez IDs are the associated names)
                        organism = "hsa", # supported organisms listed below
                        minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                        pvalueCutoff = 0.05, # padj cutoff value
                        verbose = TRUE)

# Extract the GSEA results
gseaKEGG_results.rna <- gseaKEGG.rna@result
gseaKEGG_results.rpf <- gseaKEGG.rpf@result

# Write GSEA results to file
View(gseaKEGG_results.rna)
write.csv(gseaKEGG_results.rna, "results/gseaHN_kegg.rna.csv", quote=F)
View(gseaKEGG_results.rpf)
write.csv(gseaKEGG_results.rpf, "results/gseaHN_kegg.rna.csv", quote=F)

# Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG.rpf, geneSetID = 'hsa03010')

# Output images for a single significant KEGG pathway
detach("package:dplyr", unload=TRUE) # unload dplyr to avoid conflicts
library(pathview)
hsa03010 <- pathview(gene.data = foldchanges.rpf,
                     pathway.id = "hsa03010",
                     species = "hsa",
                     limit = list(gene = max(abs(foldchanges.rpf)), # value gives the max/min limit for foldchanges
                                  cpd = 1))

## Running score and preranked list of GSEA results
gseaplot(gseaKEGG.rna, geneSetID = 1, by = "runningScore", title = gseaKEGG.rna$Description[1])
gseaplot(gseaKEGG.rna, geneSetID = 1, by = "preranked", title = gseaKEGG.rna$Description[1])
gseaplot(gseaKEGG.rna, geneSetID = 1, title = gseaKEGG.rna$Description[1])
enrichplot::gseaplot2(gseaKEGG.rna, geneSetID = 1:3)
enrichplot::gseaplot2(gseaKEGG.rna, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
enrichplot::gsearank(gseaKEGG.rna, 1, title = gseaKEGG.rna[1, "Description"])


# GSEA GO --------------------------------------------------------------------

gsego_all.rna <- gseGO(geneList = foldchanges.rna,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       keyType = "ENTREZID",
                       minGSSize = 100,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = TRUE)

gsego_all.rpf <- gseGO(geneList = foldchanges.rpf,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       keyType = "ENTREZID",
                       minGSSize = 100,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = TRUE)

golist.gsea <- vector(mode='list', length(catgenelists.gsea))
names(golist.gsea) <- names(catgenelists.gsea)
gosummary.gsea <- vector(mode='list', length(catgenelists.gsea))
names(gosummary.gsea) <- names(catgenelists.gsea)

# Create for loop to conduct GO analysis on categories and output results 
# to a data frame

for (j in 1:length(catgenelists.gsea)) {
  input <- catgenelists.gsea[[j]]
  golist.gsea[[j]] <- gseGO(geneList = input,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         keyType = "ENTREZID",
                         minGSSize = 100,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "none",
                         verbose = TRUE)
  gosummary.gsea[[j]] <- data.frame(clusterProfiler::simplify(golist.gsea[[j]]))
}

# Reactome pathway gene set enrichment analysis --------------------------------------------
# Initialize variable for for-loop
reaclist.gsea <- vector(mode='list', length(categoryRank))
names(reaclist.gsea) <- names(categoryRank)
reacsummary.gsea <- vector(mode='list', length(categoryRank))
names(reacsummary.gsea) <- names(categoryRank)

# Create for loop to conduct Reactome analysis on categories and output results to a list of dataframes
for (j in 1:length(categoryRank)) {
  input <- categoryRank[[j]]
  reaclist.gsea[[j]] <- ReactomePA::gsePathway(geneList = input, # NOTE: gene list is order ranked by RPF LFC
                                               organism = "human",
                                               minGSSize = 10,
                                               maxGSSize = 500,
                                               pvalueCutoff = 0.2,
                                               pAdjustMethod = "BH",
                                               verbose = TRUE)
  reacsummary.gsea[[j]] <- data.frame(reaclist.gsea[[j]])
}

# Pathway visualization 
library(ReactomePA)

# Chromatin modifying enzymes
cat_lm <- readr::read_tsv(file = "/Users/himanshudashora/Downloads/cat_lm.txt")
cat_mh <- readr::read_tsv(file = "/Users/himanshudashora/Downloads/cat_mh.txt")
pathway <- "HATs acetylate histones"
viewPathway(pathway,
            readable = T,
            foldChange = c(categoryRank[["L_M"]],categoryRank[["M_H"]]))
str_split(as.character(cat_lm[1,7]), ",")
str_split(as.character(cat_mh[8,7]), ",")
#HATs acetylate histones
#PKMTs methylate histone lysines

# Cholesterol biosynthesis by SREBP
pathway <- "Regulation of cholesterol biosynthesis by SREBP (SREBF)"
viewPathway(pathway,
            readable = T,
            foldChange = foldchanges.rpf)
str_split(as.character(cat_lm[30,7]), ",")

pathway <- "Signaling by Rho GTPases"
viewPathway(pathway,
            readable = T)
str_split(as.character(cat_lm[30,7]), ",")

# Signaling by NOTCH1
cat_lh <- readr::read_tsv(file = "/Users/himanshudashora/Downloads/cat_lh.txt")
pathway <- "Signaling by NOTCH1"
viewPathway(pathway,
            readable = T,
            foldChange = categoryRank[["L_H"]])
str_split(as.character(cat_lh[5,7]), ",")

# SUMOylation
pathway <- "SUMOylation of DNA damage response and repair proteins"
viewPathway(pathway,
            readable = T,
            foldChange = categoryRank[["L_H"]])
str_split(as.character(cat_lh[2,7]), ",")

# Cellular Senescence
pathway <- "Cellular Senescence"
viewPathway(pathway,
            readable = T,
            foldChange = categoryRank[["M_H"]])
str_split(as.character(cat_mh[23,7]), ",")

# Signaling by Hedgehog
pathway <- "Signaling by Hedgehog"
viewPathway(pathway,
            readable = T,
            foldChange = categoryRank[["M_H"]])
str_split(as.character(cat_mh[9,7]), ",")


# Deubiquitination
pathway <- "Ub-specific processing proteases"
pathway <- "Deubiquitination"
viewPathway(pathway,
            readable = T,
            foldChange = categoryRank[["M_H"]])
str_split(as.character(cat_mh[12,7]), ",")

# Signaling by WNT
pathway <- "Beta-catenin independent WNT signaling"
viewPathway(pathway,
            readable = T,
            foldChange = categoryRank[["L_M"]])
str_split(as.character(cat_lm[34,7]), ",")

# Signaling by RTK
pathway <- "Signaling by Receptor Tyrosine Kinases"
viewPathway(pathway,
            readable = T,
            foldChange = categoryRank[["L_M"]])
str_split(as.character(cat_mh[21,7]), ",")

cat_lm <- readr::read_tsv(file = "/Users/himanshudashora/Downloads/cat_lm.txt")

# Extract reactome genelists
require(ReactomePA)
x <- enrichPathway(names(geneList)[1:100])
geneInCategory(x)[["R-HSA-1655829.2"]]

# MSigDB ------------------------------------------------------------------

library(msigdb)
library(ExperimentHub)
library(GSEABase)

eh <- ExperimentHub()
query(eh, 'msigdb')
msigdb.hs <- getMsigdb(org = "hs", id = "EZID", version = "7.5")

length(msigdb.hs)

gs <- msigdb.hs[["BUFFA_HYPOXIA_METAGENE"]]
geneIds(gs)
collectionType(gs)
bcCategory(collectionType(gs))
bcSubCategory(collectionType(gs))
description(gs)
details(gs)

listCollections(msigdb.hs)
listSubCollections(msigdb.hs)
subsetCollection(msigdb.hs, 'h')

library(limma)

allg <- unique(unlist(geneIds(msigdb.hs)))
emat <- matrix(0, nrow = length(allg), ncol = 6)
rownames(emat) <- allg
colnames(emat) <- paste0('sample', 1:6)
head(emat)

hallmarks <- subsetCollection(msigdb.hs, 'h')
msigdb_ids <- geneIds(hallmarks)

fry_indices <- ids2indices(msigdb_ids, rownames(emat))
fry_indices[1:2]


# GSVA --------------------------------------------------------------------

library(GSVA)

p <- 10000 ## number of genes
n <- 30    ## number of samples
## simulate expression values from a standard Gaussian distribution
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
X[1:5, 1:5]

## sample gene set sizes
gs <- as.list(sample(10:100, size=100, replace=TRUE))
## sample gene sets
gs <- lapply(gs, function(n, p)
  paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))

gsva.es <- gsva(X, gs, verbose=FALSE)
dim(gsva.es)
gsva.es[1:5, 1:5]

gsva(rna_dds, msigdb.hs, verbose = TRUE)

romer()
roast()
fry()
camera()


# msigdbr -----------------------------------------------------------------

library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens")
head(all_gene_sets)
h_gene_sets <- msigdbr(species = "human", category = "H")
head(h_gene_sets)
cgp_gene_sets <- msigdbr(species = "human", category = "C2", subcategory = "CGP")
head(cgp_gene_sets)
gs <- cgp_gene_sets %>%
  dplyr::filter(gs_name == "BUFFA_HYPOXIA_METAGENE")

library(clusterProfiler)
library(tidyverse)
gene_ids_vector <- read_csv("results/rna_results.csv") %>%
  select(entrez, log2FoldChange) %>%
  dplyr::distinct(entrez, .keep_all = TRUE)
## feature 1: numeric vector
geneList <- gene_ids_vector$log2FoldChange
## feature 2: named vector
names(geneList) <- pull(gene_ids_vector, entrez)
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

msigdbr_t2g <- gs %>%
  dplyr::distinct(gs_name, entrez_gene) %>%
  as.data.frame()
enricher(gene = gene_ids_vector, TERM2GENE = msigdbr_t2g)

y <- GSEA(geneList, TERM2GENE = msigdbr_t2g)
head(y)
gseaplot(y, "BUFFA_HYPOXIA_METAGENE")

library(GSVA)
expr <- counts(rna_dds, normalized = T)
msigdbr_list <- split(x = gs$entrez_gene, f = gs$gs_name)
GSVA::gsva(expr, gset.idx.list = msigdbr_list, method = "ssgsea")
