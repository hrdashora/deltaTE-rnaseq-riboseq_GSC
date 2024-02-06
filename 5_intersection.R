### STEP 5: Intersect DTEGs from distinct datasets ###

## This script will read the list of DTEGs and determine the overlapping genes
## across the different data sets. Direction of changes will be preserved.

# Create tables of DTEG categories created in STEP 4 via the DTEG framework
## Exclusive DTEGs
exclusive_tab_gscnpc <- res_gscnpc[which(res_gscnpc$padj < 0.05 &
                                           res_gscnpc_ribo$padj < 0.05 &
                                           res_gscnpc_rna$padj > 0.05),]
exclusive_tab_gscnpc_pos <- exclusive_tab_gscnpc[which(exclusive_tab_gscnpc$log2FoldChange > 0),]
exclusive_tab_gscnpc_neg <- exclusive_tab_gscnpc[which(exclusive_tab_gscnpc$log2FoldChange < 0),]


exclusive_tab_hn <- res_hn[which(res_hn$padj < 0.05 &
                                   res_hn_ribo$padj < 0.05 &
                                   res_hn_rna$padj > 0.05),]
exclusive_tab_hn_pos <- exclusive_tab_hn[which(exclusive_tab_hn$log2FoldChange > 0),]
exclusive_tab_hn_neg <- exclusive_tab_hn[which(exclusive_tab_hn$log2FoldChange < 0),]

## Intensified & Buffered DTEGs
both_tab_gscnpc <- res_gscnpc[which(res_gscnpc$padj < 0.05 &
                                      res_gscnpc_ribo$padj < 0.05 &
                                      res_gscnpc_rna$padj < 0.05),]

both_tab_hn <- res_hn[which(res_hn$padj < 0.05 &
                              res_hn_ribo$padj < 0.05 &
                              res_hn_rna$padj < 0.05),]

intensified_tab_gscnpc <- both_tab_gscnpc[which(
  res_gscnpc[rownames(both_tab_gscnpc),2]*
    res_gscnpc_rna[rownames(both_tab_gscnpc),2] > 0),]

intensified_tab_gscnpc_pos <- intensified_tab_gscnpc[which(intensified_tab_gscnpc$log2FoldChange > 0),]
intensified_tab_gscnpc_neg <- intensified_tab_gscnpc[which(intensified_tab_gscnpc$log2FoldChange < 0),]

intensified_tab_hn <- both_tab_hn[which(
  res_hn[rownames(both_tab_hn),2]*
    res_hn_rna[rownames(both_tab_hn),2] > 0),]

intensified_tab_hn_pos <- intensified_tab_hn[which(intensified_tab_hn$log2FoldChange > 0),]
intensified_tab_hn_neg <- intensified_tab_hn[which(intensified_tab_hn$log2FoldChange < 0),]

buffered_tab_gscnpc <- both_tab_gscnpc[which(
  res_gscnpc[rownames(both_tab_gscnpc),2]*
    res_gscnpc_rna[rownames(both_tab_gscnpc),2] < 0),]

buffered_tab_gscnpc <- rbind(
  res_gscnpc[which(res_gscnpc$padj < 0.05 &
                     res_gscnpc_ribo$padj > 0.05 &
                     res_gscnpc_rna$padj < 0.05),],
  buffered_tab_gscnpc
)

buffered_tab_gscnpc_pos <- buffered_tab_gscnpc[which(buffered_tab_gscnpc$log2FoldChange > 0),]
buffered_tab_gscnpc_neg <- buffered_tab_gscnpc[which(buffered_tab_gscnpc$log2FoldChange < 0),]

buffered_tab_hn <- both_tab_hn[which(
  res_hn[rownames(both_tab_hn),2]*
    res_hn_rna[rownames(both_tab_hn),2] < 0),]

buffered_tab_hn <- rbind(
  res_hn[which(res_hn$padj < 0.05 &
                 res_hn_ribo$padj > 0.05 &
                 res_hn_rna$padj < 0.05),],
  buffered_tab_hn)

buffered_tab_hn_pos <- buffered_tab_hn[which(buffered_tab_hn$log2FoldChange > 0),]
buffered_tab_hn_neg <- buffered_tab_hn[which(buffered_tab_hn$log2FoldChange < 0),]

# Merge tables from different DTEG categories
tab_gscnpc <- rbind(exclusive_tab_gscnpc,
                    intensified_tab_gscnpc,
                    buffered_tab_gscnpc)
tab_gscnpc_pos <- tab_gscnpc[which(tab_gscnpc$log2FoldChange > 0), ]
tab_gscnpc_neg <- tab_gscnpc[which(tab_gscnpc$log2FoldChange < 0), ]

tab_hn <- rbind(exclusive_tab_hn,
                intensified_tab_hn,
                buffered_tab_hn)
tab_hn_pos <- tab_hn[which(tab_hn$log2FoldChange > 0), ]
tab_hn_neg <- tab_hn[which(tab_hn$log2FoldChange < 0), ]

# Visualize overlap with Venn diagram
## Intersect the DTEGs identified from the 'hypoxia versus normoxia' contrast with
## DTGs (RNA-seq only) identified from the 'GSC versus NPC' contrast
require(VennDiagram)

# Upregulated
A <- read.table(file = "results/DEGs_gscnpc.txt", header = TRUE, row.names = 1) %>%
  filter(log2FoldChange > 0) %>% rownames()
B <- rownames(tab_hn_pos)
overlapAB <- calculate.overlap(
  x = list("GSC v. NPC DEGs"= A, "Hypoxia v. Normoxia DTEGs"= B))

venn.plot <- draw.pairwise.venn(
  area1 = length(overlapAB$a1),
  area2 = length(overlapAB$a2),
  cross.area = length(overlapAB$a3),
  category = c("GSC v. NPC\nUpregulated DEGs in GSCs", "Hypoxia v. Normoxia\n Upregulated DTEGs in Hypoxic GSCs"), 
  fill = c("blue", "red"),
  lty = "blank",
  cex = 2,
  cat.cex = 1,
  cat.pos = c(180, 180),
  cat.dist = 0.05,
  cat.just = list(c(0, 1), c(1, 1))
)
grid.draw(venn.plot)

# Downregulated
C <- read.table(file = "results/DEGs_gscnpc.txt", header = TRUE, row.names = 1) %>%
  filter(log2FoldChange < 0) %>% rownames()
D <- rownames(tab_hn_neg)
overlapCD <- calculate.overlap(
  x = list("GSC v. NPC DEGs"= C, "Hypoxia v. Normoxia DTEGs"= D))

venn.plot <- draw.pairwise.venn(
  area1 = length(overlapCD$a1),
  area2 = length(overlapCD$a2),
  cross.area = length(overlapCD$a3),
  category = c("GSC v. NPC\nDownregulated DEGs in GSCs", "Hypoxia v. Normoxia\n Downregulated DTEGs in Hypoxic GSCs"), 
  fill = c("blue", "red"),
  lty = "blank",
  cex = 2,
  cat.cex = 1,
  cat.pos = c(180, 180),
  cat.dist = 0.05,
  cat.just = list(c(0, 1), c(1, 1))
)
grid.draw(venn.plot)

# Check Reactome pathways of Venn diagram overlap
require(ReactomePA)
de_pos <- overlapAB$a3
epath_pos <- enrichPathway(gene = de_pos,
                           organism = "human",
                           readable = TRUE)
dotplot(epath_pos,
        showCategory = 10,
        label_format = 60,
        font.size = 12)

de_neg <- overlapCD$a3
epath_neg <- enrichPathway(gene = de_neg,
                           organism = "human",
                           readable = TRUE)

dotplot(epath_neg,
        showCategory = 10,
        label_format = 60,
        font.size = 12)
