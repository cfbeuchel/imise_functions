master_regulator_network <- function(dat){
  
  # create nodes ----
  
  nodes <- melt(dat,measure.vars = c("gene", "metabolite"),
                variable.name = "group", 
                value.name = "label")
  
  # remove duplicates and unused column (only used for edges)
  nodes<- unique(nodes[,.(label, group)])
  nodes[,id := 1:.N]
  
  # create edges ----
  
  edges <- dat[, .(from = gene,
                   to = metabolite,
                   sign = ifelse(sign == 1, "pos", "neg"),
                   rsquared
                   )]
  
  # create network ----
  
  network <- graph_from_data_frame(
    d=edges, 
    vertices=nodes, 
    directed=F) 
  
  # modify properties ----
  
  # add colors
  colrs <- c("gene" = pal.genes[1],
             "metabolite" = pal.metabs[1]
             )
  V(network)$color <- colrs[V(network)$group]
  
  # shape of node
  node.shape <- c("gene" = "circle",
                  "metabolite" = "square"
                  )
  V(network)$shape <- node.shape[V(network)$group]
  
  # font of node
  node.font <- c("gene" = 1,
                 "metabolite" = 1)
  V(network)$font <- node.font[V(network)$group]
  
  # edge width based on r² of each association
  E(network)$width <- E(network)$rsquared * 200
  
  # edge color
  E(network)$edge.color <- alpha("grey35", .4)
  
  # edge line type - effect sign
  linetypes <- c("pos" = "longdash",
                 "neg" = "dotted")
  E(network)$edge.lty <- linetypes[E(network)$sign]
  
  # node label dist
  node.dist<- c(
    "gene"= .6,
    "metabolite" = -.6
  )
  V(network)$node.label.dist <- node.dist[V(network)$group]
  
  # layout ----
  
  my.layout <- layout_nicely(network)
  
  # main plot function
  plot.igraph(network,
              
              edge.color   = E(network)$edge.color,
              edge.lty     = E(network)$edge.lty,
              edge.curved  = F,
              
              vertex.color = adjustcolor(V(network)$color, alpha.f = 0.9),
              vertex.shape = V(network)$shape,
              vertex.size  = 4,
              vertex.frame.color  = NA,
              
              vertex.label.dist   = V(network)$node.label.dist,
              vertex.label.color  ="grey10",
              vertex.label.family ="Helvetica",
              vertex.label.font   =V(network)$font,
              vertex.label.cex    = .9,
              
              layout = my.layout
  )
  
  # legend ----
  
  # legend parameters
  lgd.cex <- .9
  
  # shape (+col)
  # Legend: Node
  legend(
    x="bottom", #-1,5
    c("Gene","Metabolite"), 
    pch=c(21, 22),
    col=NA, 
    pt.bg=colrs, 
    pt.cex=2, 
    cex=lgd.cex, 
    bty="n", 
    ncol=1)
  
  
  # edge lty
  legend(x="bottomright", 
         c(
           "Positive standardized effect estimate",
           "Negative standardized effect estimate"
         ), 
         lty=linetypes,
         col="grey35",
         lwd=1.5,
         cex=lgd.cex,
         bty="n", 
         ncol=1)
  
  
    thinkness <- E(network)$width %>% quantile(probs=c(0,.5,1))
  
  # edge size
  legend(x="bottomleft",
         title = c("Metabolite variance explained\nby gene expression (r²)"),
         c("0.5%",
           "0.7%",
           "1.5%"),
         lwd=thinkness,
         bty="n",
         ncol=1,
         cex=lgd.cex
  )
  
  
}