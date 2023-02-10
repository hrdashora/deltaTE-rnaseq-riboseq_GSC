# Create ORFik experiment (starting point if alignment is finished)
library(ORFik)
conf <- config.exper(experiment = "Hypoxic_Adaptation_GSC",
                     assembly = "Homo_sapiens_hg19",
                     type = c("Ribo-seq","RNA-seq"))
dir <- "/Volumes/LACIE (2 TB)/Jennifer_Project_mapped_data"
create.experiment(dir = dir,
                  exper = "Hypoxic_Adaptation_GSC",
                  fa = ,
                  txdb = ,
                  organism = ,
                  pairedEndBam = ,
                  rep = ,
                  condition = rep(,each = ))
create.experiment(file.path(),
                  exper = ,
                  fa = ,
                  txdb = ,
                  organism = ,
                  pairedEndBam = ,
                  rep = ,
                  condition = rep(,each = ))