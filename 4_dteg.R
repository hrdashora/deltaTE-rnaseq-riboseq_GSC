### STEP 4: Calculate deltaTE using the DTEG.R framework ###

## This script will utilize the protocol outlined in DTEG.R to calculate deltaTE
## and complete the appropriate visualizations

# Prepare input files
colData(se.rnaseq)
rna_counts_npc <- assay(subset(se.rnaseq, select = cell_type == "NPC"))
colnames(rna_counts_npc) <- colData(subset(se.rnaseq, select = cell_type == "NPC"))$file_name
rna_counts_npc

