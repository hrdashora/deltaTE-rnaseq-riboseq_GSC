### STEP 1: Generating files for differential expression analysis ###

## This script will prepare pre-processed BAM files for counting, then preserve
## the counts as SummarizedExperiment objects for use by DESeq2 subsequently.

# Load libraries and custom functions
suppressPackageStartupMessages({
  library(tidyverse)
  library(Rsamtools)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(GenomicAlignments)
  library(BiocParallel)  # Register multiple CPUs to Posit Server session
  #library(Rsubread)      # package not available on Posit Server
})

source("mutate_when.R")

# Prepare BAM files for counting ------------------------------------------
## Prepare BAM files from GSCs and NPCs using data from CLP and JSY Drives.
## Both locations are stored on Isilon and are soft linked to the Posit Server
## home directory under the link 'yulab'.

## Outside Posit Server environment, locations are accessible under
## "Volumes/yuj2lab/HRD/".

## Once code section has been run, save the output objects as "se.rnaseq" and
## "se.riboseq".

## CLP samples are all considered "Normoxic", but this needs to be validated

# Locate RNA-seq input files in CLP Drive
filenames.rnaseq.216 <-
  list.files(
    path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF216_RNA_2019/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE)
filenames.rnaseq.315 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF315_RNA_2019/BAM",
  pattern = "^S.*\\.bam$",
  full.names = TRUE)
filenames.rnaseq.318 <- c( # GSC 318 RNA-seq data is stored in two folders
  list.files(
    path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF318_RNA_2020/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE),
  list.files(
    path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RNA_BySampleName/GSC_CCF318_RNA_2021/BAM",
    pattern = "^S.*\\.bam$",
    full.names = TRUE)
  ) 
filenames.rnaseq.h1 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RNA_BySampleName/NPC_H1_RNA_2021/BAM",
  pattern = "^S.*\\.bam$",
  full.names = TRUE)
filenames.rnaseq.h9 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RNA_BySampleName/NPC_H9_RNA_2021/BAM",
  pattern = "^H9.*\\.bam$",
  full.names = TRUE)

filenames.rnaseq <- c(filenames.rnaseq.216,
                      filenames.rnaseq.315,
                      filenames.rnaseq.318,
                      filenames.rnaseq.h1,
                      filenames.rnaseq.h9,
                      list.files(
                        path = "/home/dashorh/yuj2lab/HRD/Sequencing/RNA Core/Jennifer_Project_mapped_data",
                        pattern = "^S0.*\\.sorted.bam$",
                        full.names = TRUE)) # add JSY Drive rna-seq files

file.exists(filenames.rnaseq) # confirm file names are valid

# Locate RIBO-seq input files in CLP Drive
filenames.riboseq.216 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF216_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.315 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF315_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.318 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RPF_BySampleName/GSC_CCF318_RPF_2020/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.h1 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RPF_BySampleName/NPC_H1_RPF_2021/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)
filenames.riboseq.h9 <- list.files(
  path = "/home/dashorh/yuj2lab/HRD/CLP Drive/RPF_BySampleName/NPC_H9_RPF_2021/BAM",
  pattern = "^rS.*\\.bam$",
  full.names = TRUE)

filenames.riboseq <- c(filenames.riboseq.216,
                       filenames.riboseq.315,
                       filenames.riboseq.318,
                       filenames.riboseq.h1,
                       filenames.riboseq.h9,
                       list.files(
                         path = "/home/dashorh/yuj2lab/HRD/Sequencing/RNA Core/Jennifer_Project_mapped_data",
                         pattern = "^rS0.*\\.bam$",
                         full.names = TRUE)) # add JSY Drive RIBO-seq files

file.exists(filenames.riboseq) # confirm file names are valid

# Create references to BAM files
bamfiles.rnaseq <- BamFileList(filenames.rnaseq, yieldSize = 2e6)
seqinfo(bamfiles.rnaseq[1]) # chromosome names for genomic features must match reference genome
bamfiles.riboseq <- BamFileList(filenames.riboseq, yieldSize = 2e6)
seqinfo(bamfiles.riboseq[1])

for(i in 1:length(filenames.rnaseq)) { # check BAM file headers in RNA-seq files
  scanBamHeader(BamFile(filenames.rnaseq[i]))
}

for(i in 1:length(filenames.riboseq)) { # check BAM file headers in RIBO-seq files
  scanBamHeader(BamFile(filenames.riboseq[i]))
}

# Obtain count matrices from BAM files ------------------------------------
## Count reads that overlap with exons of known genes.

# Define gene models
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by = "gene")

default <- registered() # store default register
register(MulticoreParam(workers = 40)) # use 40 cores; Unix and Mac only
registered() # check registered back-ends

## Run the following functions separately and save data in between
## as they take a long time to run

main_dir <- getwd() # getting project working directory
sub_dir <- "data" # set name of sub directory for data files
if (!file.exists(sub_dir)){
  # create a new sub directory inside the working directory
  dir.create(file.path(main_dir, sub_dir))
}

se.rnaseq <- summarizeOverlaps(features = ebg,
                               reads = bamfiles.rnaseq,
                               mode = "Union",
                               singleEnd = TRUE,
                               ignore.strand = TRUE)

save(se.rnaseq,
     file = "data/RNA-RIBO_SummarizedExperiment.RData")

# Certain BAM files from RIBO-seq have the same name, leading to 'duplicate' error
names <- names(bamfiles.riboseq)
names[1:4] <- stringr::str_replace(names[1:4], pattern = "genome.bam", replacement = "genome.sorted.bam")
names(bamfiles.riboseq) <- names

se.riboseq <- summarizeOverlaps(features = ebg, 
                                reads = bamfiles.riboseq,
                                mode = "Union",
                                singleEnd = TRUE,
                                ignore.strand = TRUE)
save(se.rnaseq,
     se.riboseq,
     file = "data/RNA-RIBO_SummarizedExperiment.RData")

## Move this .RData file into Isilon for use outside LRI Network with the following command:
# rsync -r -progress -size-only /home/dashorh/rnaseq-riboseq/data/RNA-RIBO_SummarizedExperiment.RData /mnt/isilon/w_stemcell/yuj2lab/HRD/Sequencing/'RNA Core'

for (param in rev(default)) { # restore the original register
  register(param)
}

registered()

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
sampleTable <- read_delim(file = file.path(dir, "meta", "gsc_npc_metatable.txt"))
colnames(sampleTable)
sampleTable <- sampleTable %>% 
  mutate(across(.cols = c(cell_line, treatment), .fns = as_factor))

sampleTable.rna <- dplyr::filter(sampleTable, assay == "rna") %>%
  dplyr::select(-c("assay", "sample"))
sampleTable.rpf <- dplyr::filter(sampleTable, assay == "ribo") %>%
  dplyr::select(-c("assay", "sample"))

# Utilize helper function to edit file names for specific rows
# sampleTable.rpf %>% mutate_when(
#   source == "CLP", list(file_name = str_replace(file_name,
#                                                 pattern = "genome.bam",
#                                                 replacement = "genome.sorted.bam")))

# Alternatively, utilize case_when within a call to the mutate function in dplyr
sampleTable.rpf <- sampleTable.rpf %>% mutate(
  file_name = case_when(
    source == "CLP" ~ str_replace(file_name,
                                  pattern = "genome.bam",
                                  replacement = "genome.sorted.bam"),
    .default = file_name))

colData(se.rnaseq) # information about the samples or experiments, currently empty
# confirm that columns of SE are in the same order as the rows of sampleTable
colnames(se.rnaseq) == sampleTable.rna$file_name
(colData(se.rnaseq) <- DataFrame(sampleTable.rna))
colData(se.rnaseq)

colData(se.riboseq)
colnames(se.riboseq) == sampleTable.rpf$file_name
(colData(se.riboseq) <- DataFrame(sampleTable.rpf))
colData(se.riboseq)

save(se.rnaseq, se.riboseq, sampleTable.rna, sampleTable.rpf, file = "data/RNA-RIBO_SummarizedExperiment.RData")

## Move this .RData file from RStudio Server into Isilon with the following command:
# rsync -r -progress -size-only /home/dashorh/rnaseq-riboseq/data/RNA-RIBO_SummarizedExperiment.RData /mnt/isilon/w_stemcell/yuj2lab/HRD/Sequencing/'RNA Core'
