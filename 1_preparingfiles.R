### STEP 1: Generating files for differential expression analysis ###

# Load libraries
library(tidyverse)
library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicAlignments)
library(BiocParallel)
library(Rsubread)
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
# Read the metadata
## NOTE: There is a possibility that treatment groups were mislabeled at the Core
## during alignment - labels are being reversed here as a correction
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
str_remove(str_sub(colnames(se.rna), 1, 3), "0") == sampleTable.rna$sample # confirm that columns of SE are in the same order as the rows of sampleTable
(colData(se.rna) <- DataFrame(sampleTable.rna))
colData(se.rna)

colData(se.rpf)
str_remove(str_sub(colnames(se.rpf), 1, 4), "0") == sampleTable.rpf$sample
(colData(se.rpf) <- DataFrame(sampleTable.rpf))
colData(se.rpf)

# # Alternative counts with featureCounts in Rsubread
# fc.rna <- featureCounts(files = filenames.rna,
#                         annot.inbuilt = "hg19",
#                         isPairedEnd = FALSE)
# fc.rpf <- featureCounts(files = filenames.rpf,
#                         annot.inbuilt = "hg19",
#                         isPairedEnd = FALSE)
# colnames(fc.rna$counts) <- sampleTable.rna$sample
# colnames(fc.rpf$counts) <- sampleTable.rpf$sample
# head(fc.rna$counts)

save(se.rna, se.rpf, sampleTable.rna, sampleTable.rpf, file = "RNA-RIBO_SummarizedExperiment.RData")
