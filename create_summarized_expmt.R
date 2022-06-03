# simulate an RNA-seq read counts table
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

# create gene locations
rowRanges <- GRanges(rep(c("chr1","chr2"), c(50,150)),
                         IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
                         strand = sample(c("+","-"), 200, TRUE),
                         feature_id = paste0("gene", 1:200))


