build_circos <- function(){
  
  circos.clear()
  
  # set gap between sectors
  circos.par("start.degree" = 90)
  circos.par("gap.after" = c(rep(2, 9), 10))
  
  # initialize the plot using the overall group as sector index
  circos.initialize(factors = names(gl),
                    xlim = cbind(rep(0, n_group), 
                                 group_size))
  
  #' ## Mark Metabolite Super-Class
  
  # tracks are the rings ordered from the most outer circle (1) to the most inner
  # circle 1+n
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.index = 1, 
    bg.border = NA, 
    track.height = 0.07,
    panel.fun = function(x, y) {
      
      # this only gives the *sector* name (in my case, AA, AC)
      nm = get.cell.meta.data("sector.index")
      
      # this returns all the regions of one sector
      r = gl[[nm]]
      
      # created lighlight for metabolite class
      highlight.sector(
        border = col.border.light,
        sector.index = nm,
        col = superpath_color[nm],
        padding = c(-.10,0,-1,0),
        text.col = col.text,
        text = nm,
        text.vjust = -0.7,
        cex=0.7,
        niceFacing = TRUE, 
        facing="bending.inside"
      )
    })
  
  #' ## Build track for No. of significant Associations
  
  # initialize plot and set max for y 
  circos.trackPlotRegion(
    ylim = c(0, log10(4000)),
    track.index = 2,
    track.height = .15,
    track.margin = c(0,0.24), # make more room for the labels and the most outer track
    cell.padding = c(0.02,1,0,1), # add a little room to the top of the y axis
    bg.col = col.bg.dark,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      # get correct indices
      nm = get.cell.meta.data("sector.index")
      r = gl[[nm]]
      n <- length(r)
      
      # get the num of sigs
      to.plot <- plot.info[match((r), metabolite),log10(fdr5hier)]
      to.plot.eta <- plot.info[match((r), metabolite), eta1]
      to.plot[is.infinite(to.plot)] <- 0
      
      
      # build coordinate system
      circos.axis(h = "bottom",
                  labels.niceFacing = TRUE,
                  major.tick = TRUE,
                  minor.ticks = 0,
                  labels = NULL,
                  major.at = 1:n,
                  direction = "inside",
                  col = col.border.dark,
                  labels.col = col.text
      )
      
      label.cols <- mclass_color[gd.class[r]]
      
      # another axis at the top with the labels
      circos.axis(h = "top",
                  sector.index = nm, 
                  major.tick = TRUE,
                  minor.ticks = 0,
                  labels = FALSE,
                  major.at = 1:n,
                  labels.cex = .5,
                  labels.facing = "outside",
                  col = col.border.dark,
                  labels.col= col.text
      )
      
      # line indicator for axis
      circos.lines(
        x = CELL_META$xlim,
        y = c(0,0),
        lwd = ".6",
        lty = "solid", 
        col = col.border.light)
      
      circos.text(
        sector.index = nm,
        x = 1:n,
        y = CELL_META$cell.ylim[2]+.4,
        labels = r,
        facing = "reverse.clockwise",
        adj = c(1,0.5),
        niceFacing = TRUE,
        col = label.cols,
        cex = .5
      )
      
      # get color for each eta
      plot.eta.cols <- col.fun.eta(to.plot.eta)
      
      # get slightly different point sizes  
      plot.eta.cex <- 1
      # plot.eta.cex <- ifelse(is.na(to.plot.eta), .4,
      #                        ifelse(to.plot <= 0.25, .4, 
      #                               ifelse(to.plot <= 0.75, .5, .7       
      #                               )))
      
      # set points for significant assocs
      circos.points(x = 1:n,
                    y = to.plot,
                    sector.index = nm,
                    pch = 19, # use different symbol for the distinction of the legend
                    cex=plot.eta.cex,
                    col =  mclass_color[gd.class[r]]
      )
      
      # draw line from origin to each point (vertical)
      for(i in 1:n){
        circos.lines(
          x = c(i,i),
          y = c(0,to.plot[i]),
          lwd = 1.5,
          lty = "solid", 
          col = mclass_color[gd.class[r[i]]]
        )
      }
      
      # draw helper lines for the axis
      for(i in c(log10(10),log10(500), log10(2000))){
        circos.lines(
          x = CELL_META$xlim,
          y = c(i,i),
          lwd = ".6",
          lty = "dashed", 
          col = col.border.light)
      }
    })
  
  #' ### Build single y-axis
  
  # build an Y-Axis
  circos.yaxis(side = "left", 
               track.index = 2, 
               sector.index = names(gl[1]),
               labels.cex = .5,
               col = col.text,
               at = c(0, log10(10), log10(100), log10(500), log10(2000)), 
               labels = c("0", "10", "0.1k", "0.5k", "2k"), 
               labels.col = col.border.dark)
  
  # Axis Label
  circos.text(x = -0.5,
              track.index = 2, 
              sector.index = names(gl[1]),
              y = CELL_META$yrange+2,
              niceFacing = TRUE, 
              facing = "reverse.clockwise",
              col = col.text, 
              cex=0.5,
              # labels = expression("log"[10]*("# Sig. Genes")))
              labels = "# Sig. Genes"
              )
  
  #' ## Build track for eta1
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.index = 3,
    track.height = 0.03,
    track.margin = c(0.005, 0.005), # circos.par("track.margin")
    cell.padding = c(0,0,0.015,0),
    bg.border = NA,
    bg.col = col.bg.light,
    panel.fun = function(x, y) {
      
      # get correct indices
      nm = get.cell.meta.data("sector.index")
      r = gl[[nm]]
      n <- length(r)
      
      # get the Info on I2
      to.plot.eta <- plot.info[match((r), metabolite),eta1]
      
      # First draw a coordiante system
      circos.axis(h = CELL_META$yrange/2, # axis in the middle of the track
                  labels.niceFacing = FALSE,
                  major.tick = TRUE,
                  minor.ticks = FALSE,
                  labels = FALSE,
                  major.at = 1:n,
                  labels.cex = .5,
                  direction = "inside",
                  # labels.facing = "reverse.clockwise",
                  col = col.border.dark
                  # labels.col = col.text
      )
      
      # get color for each eta
      # plot.eta.cols <- col.fun.eta(to.plot.eta)
      
      # get discrete color
      plot.eta.cols <- ifelse(is.na(to.plot.eta),"grey",
                              ifelse(to.plot.eta<=0.05, col.eta[1],
                                     ifelse(to.plot.eta<=0.1, col.eta[2], 
                                            ifelse(to.plot.eta<=0.2, col.eta[3], col.eta[4]))))
      
      
      # get slightly different point sizes
      plot.eta.cex <- 1
      # plot.eta.cex <- ifelse(to.plot.eta==0, .2,
      #                        ifelse(to.plot.eta <= 0.05, .4,
      #                               ifelse(to.plot.eta <= 0.1, .5,
      #                                      ifelse(to.plot.eta <= 0.2, .6, .7
      #                                      ))))
      
      # Plot a point corresponding to the eta 1
      circos.points(x = 1:n,
                    y = rep(CELL_META$yrange/2,n),
                    col = plot.eta.cols,
                    pch = 19, # solid point
                    cex = plot.eta.cex #mod
      )
      
    })
  
  # Axis Label
  circos.text(x = -0.2,
              y = CELL_META$yrange/2,
              track.index = 3, 
              sector.index = names(gl[1]),
              niceFacing = TRUE, 
              cex = 0.5,
              col = col.text, 
              labels = expression(eta[1])
  )
  
  #' ## Single Line indicating Heterogeneity of associations
  
  #' Create track that builds sector plots for each group outlining eta0 and
  #' heterogeneity and fdr sig assocs
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.index = 4,
    track.height = 0.03,
    track.margin = c(0.005, 0.005), # circos.par("track.margin")
    cell.padding = c(0,0,0.015,0),
    bg.border = NA,
    bg.col = col.bg.light,
    panel.fun = function(x, y) {
      
      # get correct indices
      nm = get.cell.meta.data("sector.index")
      r = gl[[nm]]
      n <- length(r)
      
      # get the Info on I2
      to.plot <- plot.info[match((r), metabolite),i2.q99]
      
      # First draw a coordiante system
      circos.axis(h = CELL_META$yrange/2, # axis in the middle of the track
                  labels.niceFacing = FALSE,
                  major.tick = TRUE,
                  minor.ticks = FALSE,
                  labels = FALSE,
                  major.at = 1:n,
                  labels.cex = .5,
                  direction = "inside",
                  col = col.border.dark
      )
      
      # get the color for each heterogeneity
      # plot.cols <- sapply(to.plot,function(x){
      #   res <- ifelse(is.na(x),"grey",col.fun.het(x))
      #   return(res)
      # })
      
      # discrete plot colors
      plot.cols <- ifelse(is.na(to.plot),"grey",
                          ifelse(to.plot<=0.25, col.het[1],
                                 ifelse(to.plot<=0.75, col.het[2], col.het[4])))
      
      # get slightly different 
      plot.cex <- ifelse(is.na(to.plot), .4,1)
      # plot.cex <- ifelse(is.na(to.plot), .4,
      #                    ifelse(to.plot <= 0.25, .5, 
      #                           ifelse(to.plot <= 0.75, .6, .7       
      #                           )))
      
      # Plot the 99th percentile of I2
      circos.points(x = 1:n,
                    y = rep(CELL_META$yrange/2,n),
                    col = plot.cols,
                    pch = 19, # solid point
                    cex = plot.cex #mod
      )
      
    })
  
  # Axis Label
  circos.text(x = -0.2,
              y = 0.5,
              track.index = 4, 
              sector.index = names(gl[1]),
              niceFacing = TRUE, 
              # facing = "reverse.clockwise",
              cex = 0.5,
              col = col.text, labels = expression("I"^2))

  #' ## Build track for Zero-Inflation
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.index = 5,
    track.height = 0.03,
    track.margin = c(0.005, 0.005), # circos.par("track.margin")
    cell.padding = c(0,0,0.015,0),
    bg.border = NA,
    bg.col = col.bg.light,
    panel.fun = function(x, y) {
      
      # get correct indices
      nm = get.cell.meta.data("sector.index")
      r = gl[[nm]]
      n <- length(r)
      
      # get the Info on I2
      to.plot.zinf <- plot.info[match((r), metabolite),zero.inf]
      
      # First draw a coordiante system
      circos.axis(h = CELL_META$yrange/2, # axis in the middle of the track
                  labels.niceFacing = FALSE,
                  major.tick = TRUE,
                  minor.ticks = FALSE,
                  labels = FALSE,
                  major.at = 1:n,
                  labels.cex = .5,
                  direction = "inside",
                  col = col.border.dark
      )
      
      # get color for each eta
      plot.zinf.cols <- ifelse(to.plot.zinf==0, col.zinf[1],
                               ifelse(to.plot.zinf <= 0.3, col.zinf[2],
                                      ifelse(to.plot.zinf <= 0.5, col.zinf[3], 
                                             col.zinf[4]
                                      )))
      # plot.zinf.cols <- col.fun.zinf(to.plot.zinf)
      
      # get slightly different point sizes
      plot.zinf.cex <- 1
      # plot.zinf.cex <- ifelse(to.plot.zinf==0, .2,
      #                         ifelse(to.plot.zinf <= 0.3, .3,
      #                                ifelse(to.plot.zinf <= 0.5, .5,
      #                                       ifelse(to.plot.zinf <= 0.8, .6, .7
      #                                       ))))
      
      # Plot a point corresponding to the eta 1
      circos.points(x = 1:n,
                    y = rep(0.5,n),
                    col = plot.zinf.cols,
                    pch = 19, # solid point
                    cex = plot.zinf.cex #mod
      )
    })
  
  # Axis Label
  circos.text(x = -0.5,
              y = 0.5,
              track.index = 5, 
              sector.index = names(gl[1]),
              niceFacing = TRUE, 
              cex = 0.5,
              col = col.text, 
              labels = "%Zero-\nInfl.")
  
  #' ## Mark Metabolite Class
  
  # tracks are the rings ordered from the most outer circle (1) to the most inner
  # circle 1+n
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.index = 6, 
    bg.border = NA, 
    track.margin = c(0,0.01),
    cell.padding = c(0,1,0,1),
    track.height = 0.03,
    panel.fun = function(x, y) {
      
      # this only gives the *sector* name (in my case, AA, AC)
      nm = get.cell.meta.data("sector.index")
      
      # this returns all the regions of one sector
      r = gl[[nm]]
      n <- length(r)
      
      # get colors for each group
      plot.col.mclass <- mclass_color[gd.class[r]]
      
      # build
      circos.rect(xleft = c(1:n-.3),
                  ybottom = rep(0.1,n),
                  xright = c(1:n+.3),
                  ytop = rep(0.8,n), 
                  border=col.border.light,
                  col=plot.col.mclass)
    })
  
  #' ## Links
  
  #' This code extracts and calculates positions for the links between the metabolites
  
  v_i = NULL
  v_j = NULL
  v_g1 = NULL
  v_g2 = NULL
  v_k1 = NULL
  v_k2 = NULL
  v = NULL
  for(i in 1:(n-1)) {
    for(j in seq(i+1, n)) {
      g1 = gd[rn[i]] 
      g2 = gd[rn[j]]
      r1 = gd[gd == g1]
      k1 = which(names(r1) == rn[i]) #- 0.5
      r2 = gd[gd == g2]
      k2 = which(names(r2) == rn[j]) #- 0.5
      
      v_i = c(v_i, i)
      v_j = c(v_j, j)
      v_g1 = c(v_g1, g1)
      v_g2 = c(v_g2, g2)
      v_k1 = c(v_k1, k1)
      v_k2 = c(v_k2, k2)
      v = c(v, mat[i, j])
    }
  }
  df = data.frame(i = v_i, j = v_j, g1 = v_g1, g2 = v_g2, k1 = v_k1, k2 = v_k2, v = v)
  df = df[order(abs(df$v)), ]
  
  df$lty <- ifelse(df$g1 == df$g2, "solid", "dashed")
  
  # draw the sectors
  for(i in seq_len(nrow(df))) {
    circos.link(
      sector.index1 = df$g1[i], 
      point1 = df$k1[i], 
      sector.index2 = df$g2[i], 
      point2 = df$k2[i], 
      col = col_fun(df$v[i]),
      lty = df$lty[i]
    )
  }
  
  circos.info()
  
  #### Legends
  
  # add Legend for correlation links
  # https://jokergoo.github.io/blog/html/add_legend_to_circlize.html
  lgd_pearson <- Legend(at = c(-1, -0.5, 0, 0.5, 1), 
                        col_fun = col_fun, 
                        title_position = "topleft", 
                        title = "Pearson's correlation\nr>0.9 (FDR=5%)\nof effect estimates", 
                        direction = "horizontal", 
                        # grid_width = 0.1,
                        legend_width = unit(1, "in"),
                        labels_gp = gpar(col=col.text, cex=0.9),
                        title_gp =  gpar(col=col.text.dark))
  
  # add Legend for i2
  lgd_i2 <- Legend(at = c("No signficant associations", "0 – 0.25", "0.25 – 0.75", ">0.75"),
                   type = "points", 
                   legend_gp = gpar(col= c("grey", col.het[c(1,2,4)]),
                                    cex= c(0.4, 0.4, 0.6, 0.7)),
                   title_position = "topleft", 
                   title = expression(atop("99th percentile "*I^2," of sig. associations")), 
                   direction = "horizontal",
                   # grid_width = 0.1,
                   legend_width = unit(1, "in"),
                   labels_gp = gpar(col=col.text, cex=0.9),
                   title_gp =  gpar(col=col.text.dark))
  
  # add Legend for eta1
  
  lgd_eta1 <- Legend(
    at = c("0 – 0.05", "0.05 – 0.1", "0.1 – 0.2", "> 0.2"), 
    type = "points", 
    legend_gp = gpar(col=col.eta,
                     cex=c(0.2,0.4,0.55,0.7)),
    # pch=15, # square as symbol
    title_position = "topleft", 
    title = expression(eta[1]), 
    direction = "horizontal",
    # grid_width = 0.1, 
    legend_width = unit(1, "in"),
    labels_gp = gpar(col=col.text, 
                     cex=0.9),
    title_gp =  gpar(col=col.text.dark))
  
  # add Legend for Zero Inflation
  
  lgd_zinf <- Legend(
    at = c("0%", expression(""<=30*"%"), expression(""<=50*"%"), "> 50%"),
    type = "points",
    legend_gp = gpar(col=col.zinf),
    title_position = "topleft", 
    title = "Weighted % of exess zeros\nin metabolite data", 
    direction = "horizontal",
    # grid_width = 0.1, 
    legend_width = unit(1, "in"),
    labels_gp = gpar(col=col.text, 
                     cex=0.9),
    title_gp =  gpar(col=col.text.dark)
    )
  
  # Pack the legend in two packs
  legend.left.up = packLegend(lgd_i2,lgd_pearson, 
                           row_gap = unit(1, "cm"))
  legend.left.down = packLegend(lgd_eta1, lgd_zinf,
                           row_gap = unit(1, "cm"))

  #' Draw the legends
  
  # draw the lower two legends
  pushViewport(viewport(x = 0.04,
                        y = 0.05,
                        width = 0.1,
                        height = 0.1
  )
  )
  draw(legend.left.down, 
       x = unit(1, "cm"), 
       y = unit(1, "cm"), 
       just = c("left", "bottom"))
  popViewport()
  
  # draw the upper two legends
  pushViewport(viewport(x = 0.04,
                        y = 0.76,
                        width = 0.1,
                        height = 0.1
                        )
               )
  draw(legend.left.up, 
       x = unit(1, "cm"), 
       y = unit(1, "cm"), 
       just = c("left", "bottom"))
  popViewport()
  
  # add legend for aa/ac/mix
  
  lgd_mclass <- Legend(at = c("Amino Acids", 
                              "Mix",
                              "Acylcarnitines"
                              ), 
                       type = "grid", 
                       border = col.border.light,
                       legend_gp = gpar(fill = mclass_color), 
                       title_position = "topleft", 
                       title = "Metabolite Class", 
                       nrow = 4, 
                       gap = 0.2,
                       direction = "vertical",
                       labels_gp = gpar(col=col.text,
                                        cex=0.9),
                       title_gp =  gpar(col=col.text.dark))
  
  pushViewport(
    viewport(x = 0.83, 
             y = 0.04, 
             width = 0.1, 
             height = 0.1, 
             just = c("left", "bottom")))
  grid.draw(lgd_mclass)
  upViewport()
  
}

