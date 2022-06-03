# load package and workflow template
library(systemPipeRdata)
getWorkenvir(workflow = "rnaseq",
             mydirname = "rnaseq-riboseq")
setwd("rnaseq-riboseq")

# create the workflow
sal <- SPRproject()
sal <- importWF(sal, file_path = "systemPipeRNAseq.Rmd", verbose = FALSE)

# run workflow
sal <- runWF(sal)

# workflow visualization
plotWF(sal)

# report generation
sal <- renderReport(sal)
sal <- renderLogs(sal)