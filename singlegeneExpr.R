singlegeneExpr <- function(gene_id, ddsMat.rna, ddsMat.rpf, symbol) {
  require(ggplot2)
  require(DESeq2)
  
  # Find the Entrez ID
  #gene_id <- hg19annot[hg19annot$SYMBOL == symbol, "ENTREZID"]
  
  # Save plotcounts to a data frame object
  counts.rna <- plotCounts(ddsMat.rna, gene=gene_id, intgroup=c("treatment","cell_line"), returnData=TRUE)
  counts.rpf <- plotCounts(ddsMat.rpf, gene=gene_id, intgroup=c("treatment","cell_line"), returnData=TRUE)
  
  # Plot the normalized counts
  plot.rna <- ggplot(counts.rna, aes(x = treatment, y = count)) + 
    geom_point(aes(color = cell_line), position=position_jitter(w = 0.1,h = 0)) +
    #geom_text_repel(aes(label = rownames())) + 
    theme_bw() +
    labs(x = "Treatment", y = "Normalized Counts", color = "Cell Line") +
    ggtitle(paste(symbol,"RNA-seq")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot.rpf <- ggplot(counts.rpf, aes(x = treatment, y = count)) + 
    geom_point(aes(color = cell_line), position=position_jitter(w = 0.1,h = 0)) +
    #geom_text_repel(aes(label = rownames())) + 
    theme_bw() +
    labs(x = "Treatment", y = "Normalized Counts", color = "Cell Line") +
    ggtitle(paste(symbol, "RIBO-seq")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list(rna = plot.rna, rpf = plot.rpf))
}