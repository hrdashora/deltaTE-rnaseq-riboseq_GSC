#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Ribo-seq HEK293 (2020) Investigative analysis of quality of new Ribo-seq data
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Article: https://f1000research.com/articles/9-174/v2#ref-5
# Design: Wild type (WT) vs codon optimized (CO) (gene F9)
library(ORFik)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Config
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Specify paths wanted for NGS data, genome, annotation and STAR index
# If you use local files, make a conf variable with existing directories
# Seperate Ribo-seq and RNA-seq into separate folders with type argument
conf <- config.exper(experiment = "Alexaki_Human",
                     assembly = "Homo_sapiens_GRCh38_101",
                     type = c("Ribo-seq", "RNA-seq"))
# Will create default config paths, if you want more control of where the
# data is stored, check out function config() function

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Download fastq files for experiment and rename (Skip if you have the files already)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# SRA Meta data download (work for ERA, DRA and GEO too)
study <- download.SRA.metadata("PRJNA591214", auto.detect = TRUE)
# Auto detection worked, all Ribo-seq and RNA-seq samples detected
# NOTE: Could not detect condition CO, only wild type (WT)

# Split study into (Ribo-seq / RNA-seq)
study.rfp <- study[LIBRARYTYPE == "RFP",]
study.rna <- study[LIBRARYTYPE == "RNA",]
# Download fastq files (uses SRR numbers (RUN column) from study))
# The sample_title column had good names to rename files:
download.SRA(study.rfp, conf["fastq Ribo-seq"],
             rename = study.rfp$sample_title, subset = 2000000)
download.SRA(study.rna, conf["fastq RNA-seq"],
             rename = study.rna$sample_title, subset = 2000000)

# Which organism is this, scientific name, like "Homo sapiens" or "Danio rerio"
organism <- study$ScientificName[1] # Usually you find organism here, else set it yourself
paired.end.rfp <- study.rfp$LibraryLayout == "PAIRED"
paired.end.rna <- study.rna$LibraryLayout == "PAIRED"
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Annotation (Download genome, transcript annotation and contaminants)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

annotation <- getGenomeAndAnnotation(
  organism = organism,
  genome = TRUE, GTF = TRUE,
  phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
  output.dir = conf["ref"],
  assembly_type = "primary_assembly"
)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# STAR index (index the genome and contaminants seperatly)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Remove max.ram = 20 and SAsparse = 2, if you have more than 64GB ram
index <- STAR.index(annotation, wait = TRUE, max.ram = 20, SAsparse = 2)

# Show all annotations you have made with ORFik so far 
list.genomes()
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Alignment (with depletion of phix, rRNA, ncRNA and tRNAs) & (with MultiQC of final STAR alignment)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

STAR.align.folder(conf["fastq Ribo-seq"], conf["bam Ribo-seq"], index,
                  paired.end = paired.end.rfp,
                  steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                  adapter.sequence = "auto", # Adapters are auto detected
                  trim.front = 0, min.length = 20)

STAR.align.folder(conf["fastq RNA-seq"], conf["bam RNA-seq"], index,
                  paired.end = paired.end.rna,
                  steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                  adapter.sequence = "auto", # Adapters are auto detected
                  trim.front = 0, min.length = 20)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create experiment (Starting point if alignment is finished)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFik)
create.experiment(file.path(conf["bam Ribo-seq"], "aligned/"),
                  exper = conf["exp Ribo-seq"],
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  pairedEndBam = paired.end.rfp,
                  rep = c(1,2,3,1,2,3),
                  condition = rep(c("CO", "WT"), each = 3))
create.experiment(file.path(conf["bam RNA-seq"], "aligned/"),
                  exper = conf["exp RNA-seq"],
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  pairedEndBam = paired.end.rna,
                  rep = c(1,2,3,1,2,3),
                  condition = rep(c("CO", "WT"), each = 3))

library(ORFik)
# Show the experiments you have made with ORFik so far 
list.experiments()
df.rfp <- read.experiment("Alexaki_Human_Ribo-seq")
df.rna <- read.experiment("Alexaki_Human_RNA-seq")


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Convert files and run Annotation vs alignment QC
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# General QC
ORFikQC(df.rfp)
ORFikQC(df.rna)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# P-shifting of Ribo-seq reads:
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# From ORFikQC it looks like 20, 21, 27 and 28 are candidates for Ribosomal footprints
shiftFootprintsByExperiment(df.rfp, accepted.lengths = c(20:21, 27:28))
# Now check if you are happy with shifts, these libraries have some interesting
# periodicity for read length 20 and 27, 
# it has identical amount of reads in frame 0 and 1, not optimal for ORF detection.
shiftPlots(df.rfp, output = "auto", downstream = 30) # Barplots, better details
shiftPlots(df.rfp, output = "auto", downstream = 30, type = "heatmap") # Heatmaps, better overview

# Ribo-seq specific QC
remove.experiments(df.rfp) # Remove loaded data (it is not pshifted)
RiboQC.plot(df.rfp, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
remove.experiments(df.rfp)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create heatmaps (Ribo-seq)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Pre-pshifting
heatMapRegion(df.rfp, region = c("TIS", "TTS"), shifting = "5prime", type = "ofst",
              outdir = file.path(QCfolder(df), "heatmaps/pre-pshift/"))
remove.experiments(df.rfp)
# After pshifting
heatMapRegion(df.rfp, region = c("TIS", "TTS"), shifting = "5prime", type = "pshifted",
              outdir = file.path(QCfolder(df), "heatmaps/pshifted/"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Count table analysis: TE tables
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Shifting looks good, let's make count tables of pshifted libraries:
# As a note: Correlation between count tables of pshifted vs raw libs is ~ 40% usually.
countTable_regions(df.rfp, lib.type = "pshifted", rel.dir = "pshifted")

# TE per library match
countsRFP <- countTable(df.rfp, region = "cds", type = "fpkm", collapse = FALSE, count.folder = "pshifted")
countsRNA <- countTable(df.rna, region = "mrna", type = "fpkm", collapse = FALSE)
countsTE <- (countsRFP + 1) / (countsRNA + 1) # with pseudo count
# TE per condition (WT vs CO) (collapse replicates)
countsRFP <- countTable(df.rfp, region = "cds", type = "fpkm", collapse = TRUE, count.folder = "pshifted")
countsRNA <- countTable(df.rna, region = "mrna", type = "fpkm", collapse = TRUE)
countsTE <- (countsRFP + 1) / (countsRNA + 1) # with pseudo count
# TE merged all libraries
countsRFP <- countTable(df.rfp, region = "cds", type = "fpkm", collapse = "all", count.folder = "pshifted")
countsRNA <- countTable(df.rna, region = "mrna", type = "fpkm", collapse = "all")
countsTE <- (countsRFP + 1) / (countsRNA + 1) # with pseudo count

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Differential translation analysis (condition: WT vs CO)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# The design is by default chosen by this factor: The condition column in this case
design(df.rfp, multi.factor = FALSE)
# We now run, and here get 11 unique DTEG genes
res <- DTEG.analysis(df.rfp, df.rna)
# Now let's check if the CO group overexpress the F9 Gene (ENSG00000101981):
significant_genes <- res[Regulation != "No change",]
gene_names <- txNamesToGeneNames(significant_genes$id, df.rfp)
"ENSG00000101981" %in% unique(gene_names) # TRUE
# It does, good good.
# How is it regulated ?
significant_genes[which(gene_names == "ENSG00000101981"),] # By mRNA abundance
# If you downloaded the full libraries, do this to use pshifted libraries instead.
# Not a valid result for pshifted libraries using subset
res <- DTEG.analysis(df.rfp, df.rna, design = "condition", 
                     RFP_counts = countTable(df.rfp, region = "cds", type = "summarized",
                                             count.folder = "pshifted"))
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Peak detection (strong peaks in CDS)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

peaks <- findPeaksPerGene(loadRegion(df.rfp, "cds"), reads = RFP_WT_r1, type = "max")
ORFik::windowCoveragePlot(peaks, type = "cds", scoring = "transcriptNormalized")
# The gene does not have a strong max peak in WT rep1
"ENSG00000101981" %in% peaks$gene_id # FALSE

peaks_CO <- findPeaksPerGene(loadRegion(df.rfp, "cds"), reads = RFP_CO_r1, type = "max")
# The gene does not have a strong max peak in CO rep1 either
"ENSG00000101981" %in% peaks_CO$gene_id # FALSE

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Feature table (From WT rep 1)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
cds <- loadRegion(df.rfp, "cds")
cds <- ORFik:::removeMetaCols(cds) # Dont need them
cds <- cds[filterTranscripts(df.rfp)] # Filter to sane transcripts (annotation is not perfect)
dt <- computeFeatures(cds,
                      RFP = fimport(filepath(df.rfp[6,], "pshifted")),
                      RNA = fimport(filepath(df.rna[6,], "ofst")), Gtf = df.rfp,
                      grl.is.sorted = TRUE, faFile = df.rfp,
                      weight.RFP = "score", weight.RNA = "score",
                      riboStart = 21, uorfFeatures = FALSE)
# The significant DTEGs.
dt[names(cds) %in%  significant_genes$id,]
# All genes with strong 3nt periodicity of Ribo-seq
dt[ORFScores > 5,]
# Not all genes start with ATG, possible errors in annotation
table(dt$StartCodons)
# All genes with strong start codon peak
dt[startCodonCoverage > 5,]
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Gene plotting (advanced under development!)
# (Using package that extends ORFik for interactive html plots (RiboCrypt))
# Will create interactive plot for Ribo-seq and RNA-seq sample: Wild type rep 3
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# This package also available on Bioconductor since Bioc version 3.14
# BiocManager::install("RiboCrypt")
devtools::install_github("m-swirski/RiboCrypt", dependencies = TRUE)
library(RiboCrypt)
cds <- loadRegion(df.rfp, "cds")
RiboCrypt::multiOmicsPlot_list(cds[1640], cds[1640], reference_sequence = findFa(df.rfp@fafile), reads = list(fimport(filepath(df.rna[6,], "ofst")), fimport(filepath(df.rfp[6,], "pshifted"))), 
                               ylabels = c("RNA", "RFP"), withFrames = c(F, T), frames_type = "columns")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# uORF analysis (advanced under development!)
# (using the extension package to ORFik: uORFomePipe)
# Will create a mysql database + bed12 file of uORFs with color codes + plots + files with results
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
devtools::install_github("Roleren/uORFomePipe", dependencies = TRUE)
library(uORFomePipe)
find_uORFome("/media/roler/S/data/Bio_data/projects/Alexaki_uORFome/",
             df.rfp = df.rfp, df.rna = df.rna, df.cage = NULL, biomart = NULL,
             startCodons = "ATG|CTG|TTG|GTG", BPPARAM = BiocParallel::MulticoreParam(2))

grl <- getUorfsInDb()
pred <- readTable("finalPredWithProb")$prediction
cov <- readTable("startCodonCoverage")
grl[pred == 1 & rowSums(cov) > 5]