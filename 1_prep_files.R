### STEP 1: Generating files for differential expression analysis ###

## This script will prepare pre-processed BAM files for counting, then preserve
## the counts as SummarizedExperiment objects for use by DESeq2 subsequently.

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Rsamtools)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(GenomicAlignments)
  library(BiocParallel)
  library(Rsubread)
})

# Prepare BAM files for counting ------------------------------------------
## Prepare BAM files from GSCs and NPCs using data from CLP.

## Once code section has been run, save the output objects as "se.rnaseq" and
## "se.riboseq".

## Samples are all considered "Normoxic" but this needs to be validated

# Locate RNA-seq input files
filenames.rnaseq.216 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF216_RNA_2019/BAM",
  pattern = "^S.*\\.bam$",
  full.names = TRUE)
filenames.rnaseq.315 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF315_RNA_2019/BAM",
  pattern = "^S.*\\.bam$",
  full.names = TRUE)
filenames.rnaseq.318 <- c(
  list.files(
    path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF318_RNA_2020/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE),
  list.files(
    path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF318_RNA_2021/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE)
  ) # GSC 318 RNA-seq data is stored in two folders
filenames.rnaseq.h1 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/NPC_H1_RNA_2021/BAM",
  pattern = "^S.*\\.bam$",
  full.names = TRUE)
filenames.rnaseq.h9 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/NPC_H9_RNA_2021/BAM",
  pattern = "^H9.*\\.bam$",
  full.names = TRUE)
filenames.rnaseq <- c(filenames.rnaseq.216,
                      filenames.rnaseq.315,
                      filenames.rnaseq.318,
                      filenames.rnaseq.h1,
                      filenames.rnaseq.h9)
file.exists(filenames.rnaseq) # confirm file names are valid

# Locate RIBO-seq input files
filenames.riboseq.216 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF216_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.315 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF315_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filesnames.riboseq.318 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF318_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.h1 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RPF_BySampleName/NPC_H1_RPF_2021/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.h9 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RPF_BySampleName/NPC_H9_RPF_2021/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq <- c(filenames.riboseq.216,
                       filenames.riboseq.315,
                       filesnames.riboseq.318,
                       filenames.riboseq.h1,
                       filenames.riboseq.h9)
file.exists(filenames.riboseq) # confirm file names are valid

# Create references to BAM files
bamfiles.rnaseq <- BamFileList(filenames.rnaseq, yieldSize = 2e6)
seqinfo(bamfiles.rnaseq[1]) # chromosome names for genomic features must match reference genome
bamfiles.riboseq <- BamFileList(filenames.riboseq, yieldSize = 2e6)
seqinfo(bamfiles.riboseq[1])


# Obtain count matrices from BAM files ------------------------------------
## Count reads that overlap with exons of known genes.

# Define gene models
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by = "gene")

register(MulticoreParam(workers = 6)) # use 6 cores; Unix and Mac only
se.rnaseq <- summarizeOverlaps(features = ebg,
                               reads = bamfiles.rnaseq,
                               mode = "Union",
                               singleEnd = TRUE,
                               ignore.strand = TRUE)
se.riboseq <- summarizeOverlaps(features = ebg,
                                reads = bamfiles.riboseq,
                                mode = "Union",
                                singleEnd = TRUE,
                                ignore.strand = TRUE)
register(SerialParam()) # return to single core

# Access the SummarizedExperiments
se.rnaseq
se.riboseq

head(assay(se.rnaseq)) # matrix of summarized values
head(assay(se.riboseq))
colSums(assay(se.rnaseq))
colSums(assay(se.riboseq))

rowRanges(se.rnaseq) # information about the genomic ranges
str(metadata(rowRanges(se.rnaseq)))
rowRanges(se.riboseq)
str(metadata(rowRanges(se.riboseq)))


# Read and format metadata ------------------------------------------------
## Treatment groups were mislabeled at the Core during alignment. Reverse
## labels here as a correction with the fct_recode function.

dir <- getwd()
sampleTable <- read_delim(file = file.path(dir, "meta", "metatable.txt"))
colnames(sampleTable)
sampleTable <- sampleTable[, c("sample", "cell_line", "treatment")] %>% 
  mutate(across(.cols = c(cell_line, treatment), .fns = as_factor)) %>%
  mutate(treatment = fct_recode(treatment,
                                Hypoxia = "Normoxic",
                                Normoxia = "Hypoxic")) # reversing the labels
sampleTable.rna <- dplyr::slice(sampleTable, 1:8)
sampleTable.rpf <- dplyr::slice(sampleTable, 9:16)

# Prepare BAM files
filenames.rna <- list.files(path = "/Volumes/LACIE (2 TB)/Jennifer_Project_mapped_data",
                    pattern = "^S0.*\\.bam$",
                    full.names = TRUE) #names of RNA-seq input files
file.exists(filenames.rna)
filenames.rpf <- list.files(path = "/Volumes/LACIE (2 TB)/Jennifer_Project_mapped_data",
                        pattern = "^rS0.*\\.bam$",
                        full.names = TRUE) #names of RIBO-seq input files
file.exists(filenames.rpf)

bamfiles.rna <- BamFileList(filenames.rna, yieldSize = 2e6)
seqinfo(bamfiles.rna[1])
bamfiles.rpf <- BamFileList(filenames.rpf, yieldSize = 2e6)
seqinfo(bamfiles.rpf[1])

# Define gene models
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
(ebg <- exonsBy(txdb, by = "gene"))

# Obtain count matrices
register(MulticoreParam(workers = 2)) # use 2 cores; Unix and Mac only

se.rna <- summarizeOverlaps(features = ebg,
                            reads = bamfiles.rna,
                            mode = "Union",
                            singleEnd = TRUE,
                            ignore.strand = TRUE)

se.rpf <- summarizeOverlaps(features = ebg,
                            reads = bamfiles.rpf,
                            mode = "Union",
                            singleEnd = TRUE,
                            ignore.strand = TRUE)

# Access the SummarizedExperiment
se.rna
se.rpf

head(assay(se.rna)) # matrix of summarized values
head(assay(se.rpf))
colSums(assay(se.rna))
colSums(assay(se.rpf))

rowRanges(se.rna) # information about the genomic ranges
str(metadata(rowRanges(se.rna)))
rowRanges(se.rpf)
str(metadata(rowRanges(se.rpf)))

colData(se.rna) # information about the samples or experiments, currently empty
# confirm that columns of SE are in the same order as the rows of sampleTable
str_remove(str_sub(colnames(se.rna), 1, 3), "0") == sampleTable.rna$sample
(colData(se.rna) <- DataFrame(sampleTable.rna))
colData(se.rna)

colData(se.rpf)
str_remove(str_sub(colnames(se.rpf), 1, 4), "0") == sampleTable.rpf$sample
(colData(se.rpf) <- DataFrame(sampleTable.rpf))
colData(se.rpf)

save(se.rna, se.rpf, sampleTable.rna, sampleTable.rpf, file = "RNA-RIBO_SummarizedExperiment.RData")
