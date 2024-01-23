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
## Prepare BAM files from GSCs and NPCs using data from CLP and JSY Drives.

## Once code section has been run, save the output objects as "se.rnaseq" and
## "se.riboseq".

## CLP samples are all considered "Normoxic" but this needs to be validated

# Locate RNA-seq input files in CLP Drive
filenames.rnaseq.216 <-
  list.files(
    path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF216_RNA_2019/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE)
filenames.rnaseq.315 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF315_RNA_2019/BAM",
  pattern = "^S.*\\.bam$",
  full.names = TRUE)
filenames.rnaseq.318 <- c( # GSC 318 RNA-seq data is stored in two folders
  list.files(
    path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF318_RNA_2020/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE),
  list.files(
    path = "/Volumes/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF318_RNA_2021/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE)
  ) 
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
                      filenames.rnaseq.h9,
                      list.files(
                        path = "/Volumes/yuj2lab/HRD/Sequencing/RNA Core/Jennifer_Project_mapped_data",
                        pattern = "^S0.*\\.sorted.bam$",
                        full.names = TRUE)) # add JSY Drive rna-seq files

file.exists(filenames.rnaseq) # confirm file names are valid

# Locate RIBO-seq input files in CLP Drive
filenames.riboseq.216 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF216_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.315 <- list.files(
  path = "/Volumes/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF315_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.318 <- list.files(
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
                       filenames.riboseq.318,
                       filenames.riboseq.h1,
                       filenames.riboseq.h9,
                       list.files(
                         path = "/Volumes/yuj2lab/HRD/Sequencing/RNA Core/Jennifer_Project_mapped_data",
                         pattern = "^rS0.*\\.bam$",
                         full.names = TRUE)) # add JSY Drive RIBO-seq files

file.exists(filenames.riboseq) # confirm file names are valid

# Create references to BAM files
bamfiles.rnaseq <- BamFileList(filenames.rnaseq, yieldSize = 2e6)
seqinfo(bamfiles.rnaseq[1]) # chromosome names for genomic features must match reference genome
bamfiles.riboseq <- BamFileList(filenames.riboseq, yieldSize = 2e6)
seqinfo(bamfiles.riboseq[1])

for(i in 1:length(filenames.rnaseq)) { # check BAM file headers in RNA-seq files
  #open(BamFile(filenames.rnaseq[i])) 
  scanBamHeader(BamFile(filenames.rnaseq[i]))
}

for(i in 1:length(filenames.riboseq)) { # check BAM file headers in RIBO-seq files
  #open(BamFile(filenames.riboseq[i])) 
  scanBamHeader(BamFile(filenames.riboseq[i]))
}

# Obtain count matrices from BAM files ------------------------------------
## Count reads that overlap with exons of known genes.

# Define gene models
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by = "gene")

register(MulticoreParam(workers = 10)) # use 10 cores; Unix and Mac only
registered() # check registered back-ends

## Run the following functions separately and save data in between
## as they take a long time to run
se.rnaseq_CLP <- summarizeOverlaps(features = ebg,
                                   reads = bamfiles.rnaseq[1:15],
                                   mode = "Union",
                                   singleEnd = TRUE,
                                   ignore.strand = TRUE)
save(se.rnaseq_CLP, file = "data/RNA-RIBO_SummarizedExperiment.RData")
se.riboseq_CLP <- summarizeOverlaps(features = ebg,
                                    reads = bamfiles.riboseq[1:10],
                                    mode = "Union",
                                    singleEnd = TRUE,
                                    ignore.strand = TRUE)

for (param in rev(default)) { # restore the original register
  register(param)
}

registered()

# Return to single core with register(SerialParam()) if necessary

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
