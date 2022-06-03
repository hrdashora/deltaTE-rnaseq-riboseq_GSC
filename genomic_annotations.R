# load libraries
library(AnnotationHub)
library(ensembldb)
library(dplyr)

# connect to AnnotationHub
ah <- AnnotationHub()

# query AnnotationHub
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))

# extract annotations of interest
human_ens <- human_ens[["AH100643"]]

# create a gene-level dataframe
annotations_ahb <- genes(human_ens, return.type = "data.frame") %>%
  dplyr::select(gene_id_version, gene_id, gene_name, entrezid, gene_biotype) %>%
  dplyr::filter(gene_id_version %in% res_tableHN_tb$gene)

# Wait a second, we don't have one-to-one mappings!
class(annotations_ahb$entrezid)
which(map(annotations_ahb$entrezid, length) > 1)

annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()

which(is.na(annotations_ahb$gene_name)) %>% length()

which(duplicated(annotations_ahb$gene_name)) %>% length()

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)

# How many rows does annotations_ahb have?
annotations_ahb %>% nrow()

# Return only the non-duplicated genes using indices
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]

# How many rows are we left with after removing?
annotations_ahb %>% nrow()

# Determine how many of the Entrez column entries are NA
which(is.na(annotations_ahb$entrezid)) %>%  length()

# ## Create tx2gene file
# 
# # Create a transcript dataframe
# txdb <- transcripts(human_ens, return.type = "data.frame") %>%
#   dplyr::select(tx_id, gene_id)
# txdb <- txdb[grep("ENST", txdb$tx_id),]
# 
# # Create a gene-level dataframe
# genedb <- genes(human_ens, return.type = "data.frame")  %>%
#   dplyr::select(gene_id, gene_id_version, gene_name)
# 
# # Merge the two dataframes together
# tx2gene <- inner_join(txdb, genedb)
