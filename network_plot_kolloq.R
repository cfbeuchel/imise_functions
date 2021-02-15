network_plot <- function(assocResults,
                         rSquaredColumn,
                         pColumn,
                         cohort,
                         hierarchicalNetwork = T
                         ){
  
  # pColumn <- "p.hierarchical.bonferroni"
  # my.cohort = "cohort_A"
  p.min <- 0.05
  my.cohort <- cohort
  
  # change name
  setnames(assocResults, "term", "factor", skip_absent = T)
  setnames(assocResults, rSquaredColumn, "r.squared", skip_absent = T)
  
  # create data for plot
  # plotdat2 <- assocResults[get(pColumn)<=p.min &
  #                           cohort==my.cohort, .(cohort, metab, factor, r.squared)]
  # filter based on pvalue
  plotdat2 <- assocResults[
    get(pColumn)<=p.min, 
    .(cohort, metab, factor, r.squared)]
  
  # cast for use in labels # 190312 Carl edit
  all.plotdat2 <- dcast.data.table(plotdat2, formula = metab + factor ~ cohort, value.var = "r.squared") # hier weiter rsqr ~ cohort casten
  
  # re-assign for creation of labels
  plotdat2 <- all.plotdat2
  
  # melt #TAG# - here only the wanted cohort?
  # plotdat2m <- melt(plotdat2[cohort == (my.cohort), ],
  #                   id.vars = c("cohort"),
  #                   measure.vars = c("metab", "factor"))
  plotdat2m <- melt(plotdat2, 
                    id.vars = c("cohort.max"),
                    measure.vars = c("metab", "factor"))
  
  # get nodes
  knoten2 = unique(plotdat2m[, .(label = value,
                                 group = variable )])
  
  # order nodes
  setorder(knoten2, -group) # damitlegende richtigrum ist
  
  # index
  knoten2[ ,id := 1:.N]
  
  # number factors
  plotdat2[,factor_num := knoten2[match(plotdat2$factor,
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
                       factor,
                       " - ",
                       metab,
                       "<br>Max. sig. expl. Variance: <br>",
                       "LIFE Adult: ",
                       signif(`LIFE Adult`, 3),
                       "<br>",
                       "LIFE Heart: ",
                       signif(`LIFE Heart`, 3),
                       "<br>",
                       "Sorb Study: ",
                       signif(`Sorb cohort`, 3),
                       " </p>"),
        value = cohort.max # vorher r.squared
      )])) #  #  dashes = F, smooth = T, shadow = T
  
  ecken2$smooth = F
  ecken2$shadow = T
  
  # hierarchicalNetwork
  if(hierarchicalNetwork){
    
  visNetwork(
    knoten2,
    ecken2,
    height = "1000px",
    width = "1500px") %>%
    visHierarchicalLayout(sortMethod="directed") %>%
    visNodes(font= '18px arial black') %>%
    visLegend() %>%
    visOptions(highlightNearest = T) %>%
    visPhysics(enabled = FALSE) #  %>%  visEdges(arrows = 'to')
    
  } else {
    
  # Nonhierarchical plot
  visNetwork(
    knoten2,
    ecken2,
    height = "1000px",
    width = "1500px") %>%
    visLayout(randomSeed = 12, improvedLayout = T) %>%
    visNodes(font= '32px arial black', size = 12) %>%
    visLegend(stepX = 50, stepY = 50) %>%
    visOptions(highlightNearest = T) %>%
    # visPhysics(enabled = F)
    visPhysics(enabled = T, repulsion = list("nodeDistance" = 100,
                                             "centralGravity" = 0.7,
                                             "springLength" = 150,
                                             "springConstant" = 0.9,
                                             "damping" = 0.8)) #  %>%  visEdges(arrows = 'to')
  }
}
