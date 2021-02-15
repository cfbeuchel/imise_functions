network_plot <- function(assocResults,
                         rSquaredColumn,
                         node1,
                         node2
                         ){
  
  # change name
  dat <- copy(assocResults)
  setnames(dat, rSquaredColumn, "r.squared", skip_absent = T)
  setnames(dat, old = c(node1, node2), new = c("gene", "metab"), skip_absent = TRUE)
  
  # filter based on pvalue
  plotdat2 <- dat[,.(metab, gene, r.squared)]
  
  # melt #TAG# - here only the wanted cohort?
  # plotdat2m <- melt(plotdat2[cohort == (my.cohort), ],
  #                   id.vars = c("cohort"),
  #                   measure.vars = c("metab", "factor"))
  
  plotdat2m <- melt(plotdat2, measure.vars = c("metab", "gene"))
  
  # get nodes
  knoten2 = unique(plotdat2m[, .(label = value,
                                 group = variable )])
  
  # order nodes
  setorder(knoten2, -group) # damitlegende richtigrum ist
  
  # index
  knoten2[ ,id := 1:.N]
  
  # number factors
  plotdat2[,factor_num := knoten2[match(plotdat2$gene,
                                           knoten2$label), id]]
  plotdat2[,metab_num := knoten2[match(plotdat2$metab,
                                          knoten2$label),id]]
  
  # Hierarchical plot
  ecken2 = rbind(
    unique(
      plotdat2[, .(
        from = factor_num, 
        to = metab_num, 
        color = "orange", 
        title = paste0("<p>",
                       gene,
                       " - ",
                       metab,
                       "<br>sigificantly correlating: <br>",
                       signif(r.squared, 3)),
        value = r.squared # vorher r.squared
      )])) #  #  dashes = F, smooth = T, shadow = T
  
  ecken2$smooth = F
  ecken2$shadow = T
  
  # Nonhierarchical plot
  visNetwork(
    knoten2,
    ecken2,
    height = "1000px",
    width = "1500px") %>%
    visLayout(randomSeed = 12, improvedLayout = T) %>%
    visNodes(font= '20px arial black', size = 10) %>%
    visLegend() %>%
    visOptions(highlightNearest = T) %>%
    # visPhysics(enabled = F)
    visPhysics(enabled = T, repulsion = list("nodeDistance" = 1000,
                                             "centralGravity" = 0.05,
                                             "springLength" = 1500,
                                             "springConstant" = 0.99,
                                             "damping" = 1)) #  %>%  visEdges(arrows = 'to')
}
