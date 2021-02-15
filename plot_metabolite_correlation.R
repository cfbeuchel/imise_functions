plot_metabolite_correlation <- function(
  mat,
  st1
) {
  
  # preprate annotation for color labels based on class
  m1 <- match(rownames(mat), st1$rid)
  
  # super pathways matching the matrix rows/cols
  spw <- st1[m1, metab.super.path]
  
  # create a color mapping
  cls <- c("#E16A86", "#CB7F2F", "#9F9400", "#50A315", "#00AC79", "#00AAB7", 
           "#009ADE", "#A87BE4", "#DA65C3", "#E16A86")
  mapping <- data.table(
    group = unique(spw),
    col = cls
  )
  m2 <- match(spw, mapping$group)
  
  # simple heatmap of this
  ComplexHeatmap::Heatmap(
    mat, 
    col = colorRamp2(c(-1,0,1), 
                     c("#023FA5", "#E2E2E2", "#8E063B")),
    row_names_gp = gpar(
      col = mapping[m2, col],
      cex = 0.5
    ),
    column_names_gp = gpar(
      col = mapping[m2, col],
      cex = 0.5
    ),
    cluster_rows = TRUE, 
    cluster_columns = TRUE, 
    clustering_method_columns = "ward.D2", 
    clustering_method_rows = "ward.D2",
    clustering_distance_rows = "euclidean", 
    clustering_distance_columns = "euclidean")
}