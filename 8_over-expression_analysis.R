### STEP 8: Over-expression analysis per cluster###

# Load libraries
library(ReactomePA)
library(clusterProfiler)
library(topGO)
library(org.Hs.eg.db)

# Access DESeq2 results from CSV


resTested <- read_csv("results/te_categories.csv") %>%
  mutate(category_rnarpf = as.factor(category_rnarpf), category_te = as.factor(category_te))  # subset of genes that survived independent filtering in both RNA-seq and RIBO-seq

# Extract significant gene lists from RNA-seq and RIBO-seq and prepare for ORA
sig_genes.rna <- resTested %>%
  dplyr::filter(padj.rna < 0.1 & abs(log2FoldChange.rna) >= log2(2)) %>%
  select(entrez, log2FoldChange.rna)
ORAgeneList.rna <- sig_genes.rna$log2FoldChange.rna # feature 1: numeric vector
names(ORAgeneList.rna) <- dplyr::pull(sig_genes.rna, entrez) # feature 2: named vector
ORAgeneList.rna <- sort(ORAgeneList.rna, decreasing = TRUE) # feature 3: decreasing order
de.rna <- names(ORAgeneList.rna) # fold changes not needed for ORA

sig_genes.rpf <- resTested %>%
  dplyr::filter(padj.rpf < 0.1 & abs(log2FoldChange.rpf) >= log2(2)) %>%
  select(entrez, log2FoldChange.rpf)
ORAgeneList.rpf <- sig_genes.rpf$log2FoldChange.rpf
names(ORAgeneList.rpf) <- dplyr::pull(sig_genes.rpf, entrez)
ORAgeneList.rpf <- sort(ORAgeneList.rpf, decreasing = TRUE)
de.rpf <- names(ORAgeneList.rpf)

# Identify enriched GO & Reactome pathways by over-representation analysis in RNA-seq and RIBO-seq separately
ego <- enrichGO(gene = de.rna,
                universe = as.character(dplyr::pull(resTested, entrez)),
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                readable = TRUE)
head(ego)
dotplot(ego)
x <- enrichPathway(gene = de.rna, readable = TRUE)
head(x)

# Perform ORA on the genes in each cluster
categoryList <- resTested %>%
  group_by(category_rnarpf) %>%
  group_split() %>%
  set_names(levels(resTested$category_rnarpf))

# Comparing multiple gene lists
gcSample <- sapply(categoryList, simplify = FALSE, USE.NAMES = TRUE, function(df) {
  df %>%
    dplyr::pull(entrez) # convert two-column data frame into a named vector
})

ck <- compareCluster(geneClusters = gcSample, fun = "enrichPathway")
head(ck)
enrichplot::dotplot(ck, font.size = 8, showCategory = 10, label_format = 60)
cnetplot(ck)

# GO over-representation analysis with clusterProfiler -----------------------------------------------------------
# Initialize variables for for-loop
egolist <- vector(mode = 'list', length(categoryList))
names(egolist) <- names(categoryList)
egosummary <- vector(mode = 'list', length(categoryList))
names(egosummary) <- names(categoryList)

# Create for-loop to conduct GO analysis on categories and output results to a list of dataframes
for (j in 1:length(categoryList)) {
  input <- categoryList[[j]]
  gene <- as.character(input$gene_id)
  egolist[[j]] <- enrichGO(gene = gene,
                           keyType = "ENTREZID",
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           minGSSize = 10,
                           maxGSSize = 500,
                           readable = TRUE)
  egosummary[[j]] <- data.frame(clusterProfiler::simplify(egolist[[j]]))
}

# Visualize enriched GO terms as a directed acyclic graph
goplot(egolist[["TE.down"]])

subset <- sapply(egosummary, simplify = FALSE, USE.NAMES = TRUE, function(df) {
  df %>% dplyr::filter(p.adjust < 0.001) %>%
    dplyr::select(ID, Description, GeneRatio, p.adjust, geneID)
})

# Reactome pathway over-representation analysis --------------------------------------------
# Initialize variables for for-loop
ereaclist <- vector(mode = 'list', length(categoryList))
names(ereaclist) <- names(categoryList)
ereacsummary <- vector(mode = 'list', length(categoryList))
names(ereacsummary) <- names(categoryList)

# Create for-loop to conduct Reactome analysis on categories and output results to a list of dataframes
for (j in 1:length(categoryList)) {
  input <- categoryList[[j]]
  gene <- as.character(input$gene_id)
  ereaclist[[j]] <- ReactomePA::enrichPathway(gene = gene,
                                              organism = "human",
                                              pvalueCutoff = 0.05,
                                              pAdjustMethod = "BH",
                                              qvalueCutoff = 0.2,
                                              minGSSize = 10,
                                              maxGSSize = 500,
                                              readable = TRUE)
  ereacsummary[[j]] <- data.frame(ereaclist[[j]])
}


hypolist <- ReactomePA::enrichPathway(gene = as.character(sig_genes.rna$entrez),
                                      organism = "human",
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 0.2,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      readable = TRUE)
hyposummary <- data.frame(hypolist)

# Visualize overlapping categories in Reactome pathways -------------------
idx <- vector()
termvector <- vector()
padjust <- vector()
count <- vector()
desc <- vector()

for (k in 1:length(ereacsummary)) {
  idxadd <- rep(names(ereacsummary)[k], times = length(ereacsummary[[k]]$ID))
  termvectoradd <- ereacsummary[[k]]$ID # add list of Reactome paths to a vector
  padd <- ereacsummary[[k]]$p.adjust
  countadd <- ereacsummary[[k]]$Count
  descadd <- ereacsummary[[k]]$Description
  idx <- c(idx, idxadd) # concatenate vectors
  termvector <- c(termvector, termvectoradd) 
  padjust <- c(padjust, padd)
  count <- c(count, countadd)
  desc <- c(desc, descadd)
}

termmatrix <- data.frame(idx, termvector, desc, padjust, count) %>%
  mutate_at("desc", stringr::str_remove, pattern = "Homo sapiens\r: ") %>%
  group_by(idx) %>% arrange(padjust, .by_group = TRUE)
overlapterms <- termmatrix[duplicated(termmatrix$termvector) | duplicated(termmatrix$termvector, fromLast = TRUE),] %>%
  distinct()
uniqueterms <- unique(termmatrix$termvector)

upregtermmatrix <- data.frame(idx, termvector, desc, padjust, count) %>%
  filter(idx == "TE.up") %>%
  mutate_at("desc", stringr::str_remove, pattern = "Homo sapiens\r: ") %>%
  group_by(idx) %>% arrange(padjust, .by_group = TRUE)
upregoverlap <- upregtermmatrix[duplicated(upregtermmatrix$termvector)|duplicated(upregtermmatrix$termvector, fromLast=TRUE),] %>%
  distinct()

ggplot(overlapterms, aes(x=desc, y=count, fill=idx)) +
  geom_bar(stat = "identity") +
  labs(y="Count", fill="Category")+
  scale_x_discrete(limits = unique(overlapterms$desc)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "top")



# Enrichment Heatmap ------------------------------------------------------

# Create a new data-frame that has a '1' for when a gene is part of a term
# and '0' when not
annGO <- data.frame(row.names = pull(categoryList[["TE.up"]], gene_id))

for (j in 1:dim(annGO)[1]) {
  # create a matching pattern to ensure genes match exactly
  # '^GENE,' --> Match at beginning of matching string
  # ', GENE$' --> Match at end of matching string
  # 'GENE,' --> Match between first and last gene in matching string
  gene <- categoryList[["TE.up"]]$gene_id[j]
  pattern <- paste('^', gene, '|', gene, '$|', gene, sep = '')
  for (k in 1:nrow(subset[["TE.up"]])) {
    if (any(grepl(pattern, subset[["TE.up"]]$geneID[k]))) {
      annGO[j,k] <- 1
    } else {
      annGO[j,k] <- 0
    }
  }
}
colnames(annGO) <- subset[["TE.up"]][,2]

# remove terms with no overlapping genes
annGO <- annGO [, apply(annGO, 2, mean) != 0]
# remove genes with no overlapping terms
annGO <- annGO[apply(annGO, 1, mean) != 0, ]
annGO[1:5, 1:5]


# match the order of rownames in resdf.rpf with that of annGO
resAligned <- resdf.rpf[which(as.character(resdf.rpf$gene_id) %in% rownames(annGO)),]
resAligned <- column_to_rownames(resAligned, var = "gene_id")
resAligned <- resAligned[match(rownames(annGO), rownames(resAligned)),]
all(rownames(resAligned) == rownames(annGO))

source("enrichHeatmap.R")
hmap <- enrichHeatmap(results = resAligned,
                      annotations = annGO,
                      ego.summary = egosummary[["TE.up"]],
                      hg19 = hg19annot)

draw(hmap[[1]] + hmap[[2]],
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')

hmap.rna <- enrichHeatmap(results = res_aligned200.rna,
                          annotations = annGO200.rna,
                          ego.summary = ego200_cluster_summary.rna,
                          hg19 = hg19annot)

draw(hmap.rna[[1]] + hmap.rna[[2]],
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')



# Dot Plot ----------------------------------------------------------------

dotplot(ereaclist[["TE.up"]], showCategory  = 20, font.size = 8, label_format = 60)
dotplot(ereaclist[["TE.mid"]], showCategory  = 20, font.size = 8, label_format = 60)
dotplot(ereaclist[["TE.down"]], showCategory  = 20, font.size = 8, label_format = 60)

# enrichplot::dotplot(egolist[["HH"]], showCategory  = 5, font.size = 6)
# dotplot(egolist[["HL"]], showCategory  = 5, font.size = 6)
# dotplot(egolist[["LM"]], showCategory  = 5, font.size = 6)
# dotplot(egolist[["LH"]], showCategory  = 5, font.size = 6)
# dotplot(egolist[["LL"]], showCategory  = 5, font.size = 6)
# dotplot(egolist[["LM"]], showCategory  = 5, font.size = 6)
# dotplot(egolist[["MH"]], showCategory  = 5, font.size = 6)
# dotplot(egolist[["ML"]], showCategory  = 5, font.size = 6)

# convert gene ID to Symbol
edox.rna <- setReadable(ego_all.rna, 'org.Hs.eg.db', 'ENTREZID')
edox_200.rna <- setReadable(ego_200.rna, 'org.Hs.eg.db', 'ENTREZID')
edox.rpf <- setReadable(ego_all.rpf, 'org.Hs.eg.db', 'ENTREZID')
edox_200.rpf <- setReadable(ego_200.rpf, 'org.Hs.eg.db', 'ENTREZID')
edox.ll <- setReadable(egolist[["LL"]], 'org.Hs.eg.db', 'ENTREZID')
edox.hh <- setReadable(egolist[["HH"]], 'org.Hs.eg.db', 'ENTREZID')

# To color genes by log2 fold changes, we need to extract the log2 fold
# changes from our results table creating a named vector
foldchanges.rna <- sig_genes.rna$log2FoldChange
names(foldchanges.rna) <- sig_genes.rna$gene_id
foldchanges_200.rna <- sig_200.rna$log2FoldChange
names(foldchanges_200.rna) <- sig_200.rna$gene_id

foldchanges.rpf <- sig_genes.rpf$log2FoldChange
names(foldchanges.rpf) <- sig_genes.rpf$gene_id

foldchanges_lh.rpf <- cat_lh_genes$log2FoldChange.rpf
names(foldchanges_lh.rpf) <- cat_lh_genes$gene_id

# Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(edox.rna, categorySize="pvalue", foldChange = foldchanges.rna, colorEdge = TRUE)
cnetplot(edox_200.rna, categorySize="pvalue", foldChange = foldchanges_200.rna, colorEdge = TRUE)
heatplot(edox_200.rna, foldChange = foldchanges_200.rna, showCategory = 5)

# Add similarity matrix to the termsim slot of enrichment result
edox2.rna <- enrichplot::pairwise_termsim(ego_all.rna)
edox2.rpf <- enrichplot::pairwise_termsim(ego_all.rpf)
edox2.hh <- enrichplot::pairwise_termsim(egolist[["HH"]])
edox2.ll <- enrichplot::pairwise_termsim(egolist[["LL"]])
edox2_200.rna <- enrichplot::pairwise_termsim(ego_200.rna)
enrichplot::treeplot(edox2.rna)
enrichplot::treeplot(edox2_200.rna)

# Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(edox2.rna)
emapplot(edox2.rpf)
enrichplot::upsetplot(edox2.rna)
enrichplot::upsetplot(edox2.rpf)

# If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
# HN_foldchanges <- ifelse(HN_foldchanges > 2, 2, HN_foldchanges)
# HN_foldchanges <- ifelse(HN_foldchanges < -2, -2, HN_foldchanges)
# 
# cnetplot(ego, 
#          categorySize="pvalue", 
#          showCategory = 5, 
#          foldChange=HN_foldchanges, 
#          vertex.label.font=6)












