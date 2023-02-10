## How to get the CDS and UTR sequences of genes associated with colorectal cancer
# Build a gene list
library(KEGGREST)
li <- keggList("pathway", "hsa")
ptag <- names(grep("Pathways in cancer", li, value = TRUE))
ptag

tag <- gsub("path:hsa", "", ptag)

library(KEGGgraph)
dest <- tempfile()
retrieveKGML(tag, "hsa", dest, method="auto")
crids <- as.character(parseKGML2DataFrame(dest)[,1])
crgenes <- unique(translateKEGGID2GeneID(crids))
head(crgenes)

# Identify genomic coordinates
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txbygene <- transcriptsBy(txdb, "gene")[crgenes] ## subset on cancer pathway genes
map <- relist(unlist(txbygene, use.names = FALSE)$tx_id, txbygene)
map

cds <- cdsBy(txdb, "tx")
threeUTR <- threeUTRsByTranscript(txdb)
fiveUTR <- fiveUTRsByTranscript(txdb)

txid <- unlist(map, use.names = FALSE)
cds <- cds[names(cds) %in% txid]
threeUTR <- threeUTR[names(threeUTR) %in% txid]
fiveUTR <- fiveUTR[names(fiveUTR) %in% txid]


length(txid)
length(cds)
length(threeUTR)
length(fiveUTR)

cds

# Extract sequences from BSgenome
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

threeUTR_seqs <- extractTranscriptSeqs(genome, threeUTR)
fiveUTR_seqs <- extractTranscriptSeqs(genome, fiveUTR)
cds_seqs <- extractTranscriptSeqs(genome, cds)

cds_seqs

lst3 <- relist(threeUTR_seqs, PartitioningByWidth(sum(map %in% names(threeUTR))))
lst5 <- relist(fiveUTR_seqs, PartitioningByWidth(sum(map %in% names(fiveUTR))))
lstc <- relist(cds_seqs, PartitioningByWidth(sum(map %in% names(cds))))

length(map)
table(elementNROWS(map))
table(elementNROWS(lstc))
table(elementNROWS(lst3))
names(lst3)[elementNROWS(lst3) == 0L] ## genes with no 3' UTR data
table(elementNROWS(lst5))
names(lst5)[elementNROWS(lst5) == 0L] ## genes with no 5' UTR data


