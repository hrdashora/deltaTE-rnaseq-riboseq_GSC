prefilter_genes <- function(dds, group) {
  init <- nrow(dds)
  smallestGroupSize <- min(table(colData(dds)[group]))
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  final <- nrow(dds)
  print(paste("Initial number of genes: ", init))
  print(paste("Final number of genes: ", final))
  return(dds)
}