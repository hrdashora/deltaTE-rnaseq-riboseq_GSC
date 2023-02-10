# List of ranges
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
txdf <- select(edb,
               keys=keys(edb, "GENEID"),
               columns=c("GENEID","TXID"),
               keytype="GENEID")
head(txdf,20)
table(table(txdf$GENEID))

ebt <- exonsBy(edb, by="tx")

# Compare isoforms of a gene
gid <- "ENSG00000196839"
idx <- txdf$GENEID == gid
sum(idx) # how many transcripts?

txs <- txdf$TXID[idx]
txs

ebt2 <- ebt[txs]
ebt2
ebt2[[1]]

l <- length(ebt2[[1]])
ebt2[[1]][1]
ebt2[[1]][l]

# Visualizing genes with Gviz
library(Gviz)
granges2df <- function(x) { # take a GRangesList object and turn it into a data.frame
  df <- as(x, "data.frame")
  df <- df[, c("seqnames", "start", "end", "strand", "group_name")]
  colnames(df)[1] <- "chromosome"
  colnames(df)[5] <- "transcript"
  df
}

df <- granges2df(ebt2) # make a data.frame with information about the transcripts
df$gene <- gid
grt <- GeneRegionTrack(df) # make a GeneRegionTrack

range(unlist(ebt2)) # set the range of the plot

gax <- GenomeAxisTrack()
plotTracks(list(gax,grt), chromosome = "20",
           transcriptAnnotations = "transcript",
           from = 44619522 - 10000, to = 44652233 + 10000)

df$feature <- "A" # label certain features
df$feature[6] <- "B"
grt <- GeneRegionTrack(df)

plotTracks(list(gax, grt), chromosome = "20",
           transcriptAnnotation = "transcript",
           from = 44619522 - 10000, to = 44652233 + 10000,
           A = "orange", B = "dodgerblue")

# Make a toy SummarizedExperiment
library(SummarizedExperiment)
colData <- data.frame(sample = factor(1:6),
                      condition = factor(c("A","A","B","B","C","C")),
                      treated = factor(rep(0:1, 3)))
colData

# Provide rowRanges using first 10 genes in Ensembl DB
library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86
g <- genes(txdb)
g <- keepStandardChromosomes(g, pruning.mode = "coarse")
rowRanges <- g[1:10]

# Simulate "expression" measurements
exprs <- matrix(rnorm(6 * 10), ncol = 6, nrow = 10)
se <- SummarizedExperiment(assay=list("exprs"=exprs),
                           colData = colData,
                           rowRanges = rowRanges)
se
assayNames(se)

# Add data onto the rows
mcols(se)$score <- rnorm(10)
mcols(se)

# Using the ranges of a SummarizedExperiment
query <- GRanges("1", IRanges(25000, 40000))
se.sub <- se[overlapsAny(se, query), ] # pull expression data for region of interest

# Visualizing count matrix data in a Summarized Experiment
library(DESeq2)
dds <- DESeqDataSet(rse, ~ condition + treatment)
vsd <- vst(dds, blind = FALSE)

library(matrixStats)
rv <- rowVars(assay(vsd))

anno.col <- as.data.frame(colData(vsd)[,c("condition","treatment")])
anno.col

library(pheatmap)
pheatmap(assay(vsd)[head(order(rv, decreasing = TRUE), 100), ],
         annotation_col = anno.col,
         show_rownames = FALSE, show_colnames = FALSE)

# Working with DNA sequences using Biostrings
library(Biostrings)
dna <- "AAACGCG"
class(dna)
dna <- DNAString("AAACGCG") # turn char string into DNAString object
dna <- DNAStringSet(c("AAAAA","CGCG","TCG")) # vector of DNAStrings
dna
length(dna)
width(dna)
letterFrequency(dna, "C")
letterFrequency(dna, "C", as.prob = TRUE)
letterFrequency(dna, "CG", as.prob = TRUE)

dna <- DNAStringSet(c("AACTCTC", "CTCTAAA", "AAAGAG"))
matches <- vmatchPattern("CTC", dna) # look for matches of a query in a set of strings
elementNROWS(matches)

dna <- DNAString("AAACTCAAAGAGAAATTTAAA")
pd <- PDict(c("CTC", "GAG","TTT", "AAA"))
matches <- matchPDict(pd, dna)
elementNROWS(matches)

# Working with genomes
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
seqinfo(Hsapiens)
gr <- GRanges("chr1", IRanges(1e6 + c(1,101,201), width = 100))
gr
dna <- getSeq(Hsapiens, gr)
dna

# Working with a transcriptome
library(AnnotationHub)
ah <- AnnotationHub()
display(ah)
txome <- query(ah, c("ensembl", "GRCh38", "cdna.all", "release-86"))[[1]]

library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86

seqinfo(txome)

txs <- seqnames(seqinfo(txome))
head(txs)

si <- seqinfo(txome)
idx <- 5000
txp.name <- seqnames(si)[idx]
txp.name
txp.len <- seqlengths(si)[idx]
txp.len
txp.seq <- getSeq(txome, GRanges(seqnames(si)[idx], IRanges(1, txp.len)))
txp.seq

library(EnsDb.Hsapiens.v86) # use extractTranscriptSeqs function
edb <- EnsDb.Hsapiens.v86

if (file.exists("ebt.rda")) {
  load("ebt.rda")
} else {
  ebt <- exonsBy(edb, by = "tx") # exons-by-transcript
  save(ebt, file = "ebt.rda") # save the file for future use
}

txp.name2 <- sub("\\..*","", txp.name) # cut version number
txp <- ebt[txp.name2]
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
seqlevelsStyle(genome) <- "Ensembl"
txp.seq2 <- extractTranscriptSeqs(genome, txp)
txp.seq2

txp.seq == txp.seq2

