# Load a TxDb annotation package
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene # shorthand (for convenience)
seqlevels(txdb) <- seqlevels0(txdb) # reset to the original seqlevels
txdb

# Retrieving data using the select method
library(dplyr)
keys <- merge2 %>% dplyr::pull(gene_id) %>% as.character()
#keys.hl <- merge %>% dplyr::filter(category.rna.rpf == "H_L") %>% dplyr::pull(gene_id) %>% as.character()
columns(txdb)
keytypes(txdb)
tx_lens <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
txbygene <- transcriptsBy(txdb, "gene")
gene_lens <- max(width(txbygene))
tx_df <- AnnotationDbi::select(txdb,
                               keys = keys,
                               columns = c("TXNAME","TXID","TXSTRAND","TXSTART","TXEND"),
                               keytype = "GENEID") %>% filter(!is.na(TXID))

geneid <- tx_df %>% distinct(GENEID) %>% pull(GENEID)
txbygene <- txbygene[geneid] # subset on category genes
map <- relist(unlist(txbygene, use.names = FALSE)$tx_id, txbygene) 
map # IntegerList of transcript IDs for each gene in KEYS

# Obtain coordinate information for all transcripts
GR <- transcripts(txdb)
# Obtain strand information from the transcripts
tx_strand <- strand(GR)
# Compute a GRanges object that spans the promoter region
PR <- promoters(txdb, upstream = 2000, downstream = 400)
# Retrieve genomic coordinates for exons and coding sequences
EX <- exons(txdb)
CR <- cds(txdb)

# Extract quantities and group
cds <- cdsBy(txdb, by = c("tx"))
utr5 <- fiveUTRsByTranscript(txdb)
utr3 <- threeUTRsByTranscript(txdb)

txid <- unlist(map, use.names = FALSE)
cds <- cds[names(cds) %in% txid]
threeUTR <- utr3[names(utr3) %in% txid]
fiveUTR <- utr5[names(utr5) %in% txid]

length(txid)
length(cds)
length(threeUTR)
length(fiveUTR)

cds

# Convert from ranges to actual sequence
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
tx_seqs <- extractTranscriptSeqs(genome, txdb,
                                  use.names = TRUE) # RNA sequence of the spliced transcripts
utr3_seqs <- extractTranscriptSeqs(genome, threeUTR)
utr5_seqs <- extractTranscriptSeqs(genome, fiveUTR)
cds_seqs <- extractTranscriptSeqs(genome, cds)

lst3 <- relist(utr3_seqs, PartitioningByWidth(sum(map %in% names(threeUTR))))
lst5 <- relist(utr5_seqs, PartitioningByWidth(sum(map %in% names(fiveUTR))))
lstc <- relist(cds_seqs, PartitioningByWidth(sum(map %in% names(cds))))

length(map)
table(elementNROWS(map))
table(elementNROWS(lstc))
table(elementNROWS(lst3))
names(lst3)[elementNROWS(lst3) == 0L] ## genes with no 3' UTR data
table(elementNROWS(lst5))
names(lst5)[elementNROWS(lst5) == 0L] ## genes with no 5' UTR data

# Translate the coding regions
translate(cds_seqs)

# Calculate k-mer frequency in 5' UTR
library(Biostrings)
sixMers <- oligonucleotideFrequency(utr5_seqs, width = 6, with.labels = TRUE) # tabulate frequency of 6-mers in 5' UTR
fiveMers <- oligonucleotideFrequency(utr5_seqs, width = 5, with.labels = TRUE) 
fourMers <- oligonucleotideFrequency(utr5_seqs, width = 4, with.labels = TRUE)
threeMers <- trinucleotideFrequency(utr5_seqs, with.labels = TRUE)
twoMers <- dinucleotideFrequency(utr5_seqs, with.labels = TRUE)
oneMers <- oligonucleotideFrequency(utr5_seqs, width = 1, with.labels = TRUE)

# Put all data frames into list
kMers_list <- list(sixMers, fiveMers, fourMers, threeMers, twoMers, oneMers)

# Merge all data frames in list
# kMers <- Reduce(function(x,y){
#   ans <- merge(x, y, by = 'row.names', all=TRUE)
#   row.names(ans) <- ans[,"Row.names"]
#   ans[,!names(ans) %in% "Row.names"]
#   }, kMers_list) # this step takes long
kMers <- bind_cols(kMers_list)
kMers$TXID <- as.numeric(names(fiveUTR)) # add transcript id labels 



# Select the transcript IDs that is the longest for each gene
library(dplyr)
param <- left_join(tx_df, tx_lens, by = c("TXID" = "tx_id", "TXNAME" = "tx_name", "GENEID" = "gene_id")) %>%
  group_by(GENEID) %>% dplyr::slice(which.max(tx_len)) %>% mutate_at(c("TXID"), as.numeric)
items <- param %>% filter(!(utr5_len == 0)) %>% dplyr::pull(TXID) %>% as.character(.)

# G-quadruplex detection
library(pqsfinder)
pqs <- mapply(pqsfinder, utr5_seqs[items], verbose = FALSE) # high-quality PQS with defects in G-runs
lengths(pqs)

# Find ORFs in 5' UTRs
library(ORFik)
utr5_orf <- findMapORFs(fiveUTR[items], utr5_seqs[items], groupByTx = TRUE)
lengths(utr5_orf)

# Codon usage analysis of 5' UTR
library(coRdon)
ct <- coRdon::codonTable(utr5_seqs[items])
cc <- codonCounts(ct)
id <- getID(ct)
l <- getlen(ct)
head(cc)
head(l)

lengths <- as.data.frame(getlen(ct))
colnames(lengths) <- "length"
ggplot(lengths, aes(length)) + # inspect the distribution of sequence lengths
  geom_density() +
  geom_vline(xintercept = 80, colour = "red") +
  theme_light()

enc <- ENC(ct, filtering = "soft")
scuo <- SCUO(ct, filtering = "soft")

codonparams <- dplyr::bind_cols(id, enc, scuo, l)
colnames(codonparams)[1:4] <- c("TXID", "enc", "scuo", "codon.length")
codonparams <- mutate_at(codonparams, c("TXID"), as.numeric)

# gcb <- GCB(ct, seed = c(rep(TRUE, 6873), rep(FALSE, 6872))) # sequences as a seed
# halfct <- codonTable(cc[1:6873,]) # specify codonTable as a seed
# gcb2 <- GCB(ct, seed = halfct)
# all.equal(gcb, gcb2) # TRUE

# Format parameter table
library(dplyr)
library(tibble)
utr5_orflengths <- tibble::rownames_to_column(data.frame(lengths(utr5_orf)),"TXID") %>% mutate_at(c("TXID"),  as.numeric)
colnames(utr5_orflengths)[2] <- "orf.length"
utr5_pqslengths <- tibble::rownames_to_column(data.frame(lengths(pqs)),"TXID") %>% mutate_at(c("TXID"),  as.numeric)
colnames(utr5_pqslengths)[2] <- "pqs.length"
param_list <- list(param, utr5_orflengths, utr5_pqslengths, codonparams, kMers)

library(purrr)
tes <- dplyr::select(te, gene_id, treatment, log2.te, te) %>%
  pivot_wider(names_from = treatment, values_from = c(log2.te, te))
cats <- merge2 %>%
  dplyr::select(gene_id, symbol, log2.delta_te, delta_te, level.te) %>%
  left_join(., tes, by = "gene_id")



params <- param_list %>% reduce(left_join, by="TXID") %>% left_join(., cats, by = c("GENEID" = "gene_id"))
library(readr)
write_csv(params,
          file = "data/params.csv")

#save(, file="parameters.RData")


