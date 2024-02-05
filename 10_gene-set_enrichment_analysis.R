### STEP 9: Gene Set Enrichment Analysis of Reactome and MSigDB Pathways ###

# Load libraries
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(ReactomePA)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(limma)
library(edgeR)
library(RColorBrewer)

# Utilize ranked geneList from previous script to run GSEA for Reactome
head(geneList_rna)
head(geneList_rpf)

# Create ranked a geneList for all clusters
categoryRank <- sapply(categoryList, simplify = FALSE, USE.NAMES = TRUE, function(df){
  df %>%
    dplyr::select(entrez, log2FoldChange.rpf) %>% # gene lists will be ranked by RPF in the category analysis
    arrange(desc(log2FoldChange.rpf)) %>%
    deframe() # convert two-column data frame into a named vector
})

# Identify enriched Reactome pathways by GSEA in RNA-seq and RIBO-seq separately
epath_rna <- gsePathway(geneList = geneList_rna,
                        organism = "human",
                        verbose = TRUE)
head(epath_rna)

epath_rpf <- gsePathway(gene = geneList_rpf,
                        organism = "human",
                        verbose = TRUE)
head(epath_rpf)

# Initialize variables for for-loop
reaclist_gsea <- vector(mode = 'list', length(categoryRank))
names(reaclist_gsea) <- names(categoryRank)
reacsummary_gsea <- vector(mode = 'list', length(categoryRank))
names(reacsummary_gsea) <- names(categoryRank)

# Create for loop to conduct Reactome analysis on categories and output results as a list of dataframes
## NOTE: Since the ranked geneList in each cluster does not spanning the entire genome, GSEA is not as informative as ORA
for (j in 1:length(categoryRank)) {
  input <- categoryRank[[j]]
  reaclist_gsea[[j]] <- ReactomePA::gsePathway(geneList = input, # NOTE: gene list is ranked by RPF LFCs
                                               organism = "human",
                                               verbose = TRUE)
  reacsummary_gsea[[j]] <- data.frame(reaclist_gsea[[j]])
}

# Individual pathway visualization 
pathway <- "Eukaryotic Translation Elongation"
viewPathway(pathway,
            organism = "human",
            readable = T,
            foldChange = geneList_rpf)
pathway <- "Eukaryotic Translation Termination"
viewPathway(pathway,
            organism = "human",
            readable = T,
            foldChange = geneList_rpf)

# Extract Reactome pathway genelists
geneInCategory(epath_rpf)[["R-HSA-156902"]]


# MSigDB ------------------------------------------------------------------
# Choose interesting gene sets on MSigDB 
region_sets <- c("LEIN_ASTROCYTE_MARKERS",
                 "LEIN_NEURON_MARKERS",
                 "LEIN_OLIGODENDROCYTE_MARKERS")

hypoxia_sets <- c("SEMENZA_HIF1_TARGETS",
                  "BUFFA_HYPOXIA_METAGENE",
                  "WINTER_HYPOXIA_UP",
                  "WINTER_HYPOXIA_DN",
                  "WINTER_HYPOXIA_METAGENE")

glioma_sets <- c("VERHAAK_GLIOBLASTOMA_NEURAL",
                 "VERHAAK_GLIOBLASTOMA_PRONEURAL",
                 "VERHAAK_GLIOBLASTOMA_CLASSICAL",
                 "VERHAAK_GLIOBLASTOMA_MESENCHYMAL",
                 "YAMANAKA_GLIOBLASTOMA_SURVIVAL_UP",
                 "YAMANAKA_GLIOBLASTOMA_SURVIVAL_DN",
                 "BEIER_GLIOMA_STEM_CELL_UP",
                 "BEIER_GLIOMA_STEM_CELL_DN",
                 "TCGA_GLIOBLASTOMA_COPY_NUMBER_UP",
                 "TCGA_GLIOBLASTOMA_COPY_NUMBER_DN",
                 "ZHENG_GLIOBLASTOMA_PLASTICITY_UP",
                 "ZHENG_GLIOBLASTOMA_PLASTICITY_DN")

# Download data from the msigdb R package
eh <- ExperimentHub()
query(eh, 'msigdb')
msigdb.hs <- getMsigdb(org = "hs", id = "EZID", version = "7.5")

# Access and explore the data
length(msigdb.hs)
listCollections(msigdb.hs)
listSubCollections(msigdb.hs)

# Find our terms of interest in the GeneSetCollection object
gs <- msigdb.hs[c(glioma_sets, hypoxia_sets, region_sets)]
gs_ids <- geneIds(gs)

# Retrieve curated and hallmark collections
c2 <- subsetCollection(msigdb.hs, 'c2')
c2_ids <- geneIds(c2)
hallmarks <- subsetCollection(msigdb.hs, 'h')
h_ids <- geneIds(hallmarks)

# Create dummy expression data
allg <- unique(unlist(geneIds(msigdb.hs)))
emat <- matrix(0, nrow = length(allg), ncol = 6)
rownames(emat) <- allg
colnames(emat) <- paste0('sample', 1:6)
head(emat)

# Match dummy data format with rlog data from DESeq2
exprs_mat <- assay(rld.rna)
rownames(exprs_mat)
colnames(exprs_mat) <- colData(dds.rna)$sample
head(exprs_mat)

# Convert gene sets into a list of gene indices
h_indices <- ids2indices(h_ids, rownames(exprs_mat)) # index vectors for hallmark gene sets
c2_indices <- ids2indices(c2_ids, rownames(exprs_mat)) # index vectors for curated gene sets
gs_indices <- ids2indices(gs_ids, rownames(exprs_mat)) # index vectors for interesting gene sets
c2_indices[1:2]

# Create a design matrix and contrasts
cell_line <- colData(dds.rna)$cell_line
treatment <- colData(dds.rna)$treatment
design <- model.matrix(~0+treatment+cell_line)
colnames(design) <- gsub("treatment", "", colnames(design))
design
contr.matrix <- makeContrasts(
  HypovsNorm = Hypoxia - Normoxia, 
  levels = colnames(design))
contr.matrix

# Apply the camera method to the hallmark gene sets
h_cam.HypovsNorm <- camera(y = exprs_mat,
                                  index = h_indices,
                                  design = design,
                                  contrast = contr.matrix[,1])
c2_cam.HypovsNorm <- camera(y = exprs_mat,
                            index = c2_indices,
                            design = design,
                            contrast = contr.matrix[,1])
gs_cam.HypovsNorm <- camera(y = exprs_mat,
                            index = gs_indices,
                            design = design,
                            contrast = contr.matrix[,1])
head(c2_cam.HypovsNorm,5)

# Make a barcodeplot for particular gene signatures
vfit <- lmFit(exprs_mat, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix[,1])
efit <- eBayes(vfit)

barcodeplot(efit$t[,1],
            index = c2_indices$WINTER_HYPOXIA_UP, 
            index2 = c2_indices$WINTER_HYPOXIA_DN,
            main = "HypovsNorm")

barcodeplot(efit$t[,1],
            index = c2_indices$YAMANAKA_GLIOBLASTOMA_SURVIVAL_UP, 
            index2 = c2_indices$YAMANAKA_GLIOBLASTOMA_SURVIVAL_DN,
            main = "HypovsNorm")

# Use GSVA for single sample gene set enrichment, starting with gene expression data matrix
exprs_mat[1:5, 1:5] # rlog-transformed gene expression matrix
assay(dds.rna) # un-transformed gene expression from SummarizedExperiment

# Calculate GSVA enrichment scores
rna_gsva.es <- gsva(expr = dds.rna, # provide SummarizedExperiment object
                    gset.idx.list = c2, # provide GeneSetCollection object
                    annotation = "counts", # select the assay containing the input data
                    method = "gsva",
                    kcdf = "Poisson", # input expression values are integer counts (not log transformed)
                    verbose = T)

rpf_gsva.es <- gsva(expr = dds.rpf, # provide SummarizedExperiment object
                    gset.idx.list = c2, # provide GeneSetCollection object
                    annotation = "counts", # select the assay containing the input data
                    method = "gsva",
                    kcdf = "Poisson", # input expression values are integer counts (not log transformed)
                    verbose = T)

# Assess how gene expression profiles correlate between RNA-seq and RIBO-seq data
rna_lcpms <- cpm(assay(dds.rna), log = TRUE)
rpf_lcpms <- cpm(assay(dds.rpf), log = TRUE)

# Ensure gene expression profiles from both assays are same dimentions and order
rna_lcpms <- rna_lcpms[which(rownames(rna_lcpms) %in% rownames(rpf_lcpms)), ]
dim(rna_lcpms)
rpf_lcpms <- rpf_lcpms[which(rownames(rpf_lcpms) %in% rownames(rna_lcpms)), ]
dim(rpf_lcpms)

# Calculate Spearman correlations between gene expression profiles
genecorrs <- sapply(1:nrow(rna_lcpms),
                    function(i, exprnaseq, exprpfseq) cor(exprnaseq[i, ],
                                                          exprpfseq[i, ],
                                                          method = "spearman"),
                    rna_lcpms, rpf_lcpms)
names(genecorrs) <- rownames(rna_lcpms)
hist(genecorrs, xlab="Spearman correlation", main="Gene level\n(RNA-seq log-CPMs vs RIBO-seq log-CPMs)",
     xlim=c(-1, 1), col="grey", las=1)

# Calculate Spearman correlations between GSVA enrichment scores
pwycorrs <- sapply(1:nrow(assay(rna_gsva.es)),
                   function(i, esrnaseq, esrpfseq) cor(esrnaseq[i, ],
                                                       esrpfseq[i, ],
                                                       method = "spearman"),
                   assay(rna_gsva.es), assay(rpf_gsva.es))
names(pwycorrs) <- rownames(assay(rna_gsva.es))
hist(pwycorrs, xlab="Spearman correlation", main="Pathway level\n(GSVA enrichment scores)",
     xlim=c(-1, 1), col="grey", las=1)

# Heatmap of GSVA scores for msigdb signatures
treatmentOrder <- c("Hypoxia", "Normoxia")
sampleOrderByTreatment <- sort(match(rna_gsva.es$treatment, treatmentOrder),
                               index.return = TRUE)$ix
treatmentXtable <- table(rna_gsva.es$treatment)
treatmentColorLegend <-  c(Hypoxia = "red", Normoxia = "blue")
geneSetOrder <- c(glioma_sets, hypoxia_sets, region_sets)
geneSetLabels <- gsub("_", " ", geneSetOrder)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]

heatmap(assay(rna_gsva.es)[geneSetOrder, sampleOrderByTreatment],
        Rowv = NA,
        Colv = NA,
        scale = "row",
        margins = c(3,5),
        col = hmcol,
        ColSideColors = rep(treatmentColorLegend[treatmentOrder],
                            times = treatmentXtable[treatmentOrder]),
        labCol = "",
        rna_gsva.es$treatment[sampleOrderByTreatment],
        labRow = paste(toupper(substring(geneSetLabels, 1, 1)),
                       substring(geneSetLabels, 2), sep = ""),
        cexRow = 2, 
        main = " \n ")

par(xpd=TRUE)
text(0.23,1.21, "Hypoxia", col="red", cex=1.2)
text(0.36,1.21, "Normoxia", col="blue", cex=1.2)
mtext("Gene sets", side=4, line=0, cex=1.5)
mtext("Samples          ", side=1, line=4, cex=1.5)

# Conduct a differential expression analysis at a pathway level
design <- model.matrix(~0+treatment+cell_line)
colnames(design) <- gsub("treatment", "", colnames(design))
design
contr.matrix <- makeContrasts(
  HypovsNorm = Hypoxia - Normoxia, 
  levels = colnames(design))
contr.matrix

fit <- lmFit(assay(rna_gsva.es), design)
fit <- contrasts.fit(fit, contrasts = contr.matrix[,1])
fit <- eBayes(fit)
res <- decideTests(fit, p.value = 0.05) 
summary(res) # ~500 MSigDB C2 differentially expressed pathways with FDR < 5%

# Show a volcano plot of the expression changes
tt <- topTable(fit, coef = 1, n = Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
     main="", xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~Raw~P-value))
abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.05])), col=grey(0.5), lwd=1, lty=2)
points(tt$logFC[match(DEpwys, rownames(tt))],
       -log10(tt$P.Value[match(DEpwys, rownames(tt))]), pch=".", cex=5, col="darkred")
text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.05])), "5% FDR", pos=3)

# Heatmap of GSVA enrichment scores for the differentially expressed pathways -- YOU ARE HERE
DEpwys_es <- assay(rna_gsva.es)[DEpwys, ]
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("Hypoxia", "Normoxia")
sample.color.map <- colorLegend[as.character(colData(rna_gsva.es)$treatment)]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")
heatmap(DEpwys_es, ColSideColors=sample.color.map, xlab="samples",
        ylab="Pathways", margins=c(2, 20),
        labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "",
                                          rownames(DEpwys_es))), 1, 35),
        labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering))
legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")

# Use MSigDBR to perform GSEA with clusterProfiler
all_gene_sets <- msigdbr(species = "Homo sapiens")
head(all_gene_sets)
h_gene_sets <- msigdbr(species = "human", category = "H")
head(h_gene_sets)
cgp_gene_sets <- msigdbr(species = "human", category = "C2", subcategory = "CGP")
head(cgp_gene_sets)

msigdbr_t2g <- cgp_gene_sets %>%
  dplyr::distinct(gs_name, entrez_gene) %>%
  as.data.frame()

y <- GSEA(geneList_rna,
          TERM2GENE = msigdbr_t2g)
head(y)
gseaplot(y, "BUFFA_HYPOXIA_METAGENE")



