enrichHeatmap <- function(results, annotations, ego.summary, hg19) {
  require(ComplexHeatmap)
  require(circlize)
  
  # 1. Create heatmap annotations
  
  # color bar for -log10(adjusted p-value) for each gene from DE analysis
  dfMinusLog10FDRGenes <- data.frame(-log10(
    results[which(rownames(results) %in% rownames(annotations)), 'padj']))
  dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0
  
  # color bar for fold changes
  dfFoldChangeGenes <- data.frame(
    results[which(rownames(results) %in% rownames(annotations)), 'log2FoldChange'])
  
  # merge both
  dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
  colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
  dfGeneAnno[,2] <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                           ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
  colours <- list(
    'Log2FC' = c('Up-regulated' = 'royalblue', 'Down-regulated' = 'yellow'))
  
  haGenes <- rowAnnotation(
    df = dfGeneAnno,
    col = colours,
    width = unit(1,'cm'),
    annotation_name_side = 'top')
  
  # colour bar for -log10(Benjamini enrichment Q value) for GO results
  dfMinusLog10BenjaminiTerms <- data.frame(-log10(
    ego.summary[which(ego.summary$Description %in% colnames(annotations)), 'p.adjust']))
  colnames(dfMinusLog10BenjaminiTerms) <- 'Enrichment\nterm score'
  haTerms <- HeatmapAnnotation(
    Term = anno_text(
      colnames(annotations),
      rot = 45,
      just = 'right',
      gp = gpar(fontsize = 8)),
    annotation_height = unit.c(unit(8, 'cm')),
    annotation_name_side = 'left')
  haScore <- HeatmapAnnotation(
    df = dfMinusLog10BenjaminiTerms,
    simple_anno_size = unit(1, 'cm'),
    annotation_name_side = 'left')
  ha <- c(haScore, haTerms)
  
  # 2. Generate the heatmap
  idx <- match(rownames(annotations), hg19$ENTREZID)
  rownames(annotations) <- hg19$SYMBOL[idx]
  hmapGSEA <- Heatmap(as.matrix(annotations),
                      name = 'RNA GO enrichment',
                      split = dfGeneAnno[,2],
                      
                      col = c('0' = 'white', '1' = 'forestgreen'),
                      
                      rect_gp = gpar(col = 'grey85'),
                      
                      cluster_rows = TRUE,
                      show_row_dend = TRUE,
                      row_title = 'Top Genes',
                      row_title_side = 'left',
                      row_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                      row_title_rot = 90,
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 6),
                      row_names_side = 'left',
                      row_dend_width = unit(35, 'mm'),
                      
                      cluster_columns = TRUE,
                      show_column_dend = TRUE,
                      column_title = 'Enriched terms',
                      column_title_side = 'top',
                      column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                      column_title_rot = 0,
                      show_column_names = FALSE,
                      
                      show_heatmap_legend = FALSE,
                      
                      clustering_distance_columns = 'euclidean',
                      clustering_method_columns = 'ward.D2',
                      clustering_distance_rows = 'euclidean',
                      clustering_method_rows = 'ward.D2',
                      
                      bottom_annotation = ha)
  
  return(list(gsea = hmapGSEA, genes = haGenes))
}