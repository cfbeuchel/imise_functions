quotient_plot <- function(){
  
  # set plotting colors
  bg.col <- "grey97"
  axis.col <- "grey30"
  text.col <- "grey29"
  
  # set axis limits
  y.lim <- c(0, 3200)
  x.lim <- c(0,3)
  x.at <- c(.5,1.5,2.5)
  # y.at <- seq(y.lim[1], y.lim[2], by = 400) # graphics::axTicks(2)
  y.at <- c(0, 1, 10, 50, 100, 200, 400, 800, 1600, 3200)
  x.lab <- c("Amino Acid", "Quotient", "Acylcarnitine")
  y.lab <- TRUE
  
  # for log ticks
  y.lim <- c(-0.5, 3.5) # ~log10(2200)
  y.lab <- y.at
  y.at <- log10(y.at)
  y.at[1] <- -0.5
  
  plot.new()
  
  par(xpd = TRUE)
  par(mar = c( 5.1 + 5, 4.1, 4.1, 2.1))
  
  plot.window(xlim = x.lim, 
              ylim = y.lim,
              oma=c(4, 0, 0, 0))
  rect(col = "grey93", 
       xleft = par("usr")[1], 
       ybottom = par("usr")[3],
       xright = par("usr")[2],
       ytop = par("usr")[4], 
       border = NA)
  
  # draw a second rectangle to mark of the x-axis
  rect(col = bg.col, 
       xleft = 1, 
       xright = 2,
       ybottom = par("usr")[3],
       ytop = par("usr")[4], 
       border = NA)
  rect(col = bg.col, 
       xleft = -0.12, 
       xright = 0,
       ybottom = par("usr")[3],
       ytop = par("usr")[4], 
       border = NA)
  rect(col = bg.col, 
       xleft = 3, 
       xright = 3.12,
       ybottom = par("usr")[3],
       ytop = par("usr")[4], 
       border = NA)
  
  # graphics::box(bty="l",
                # lwd=2)
  
  # draw helper lines from each y tick across the plotting area
  helper.x <- c(par("usr")[1], par("usr")[2])
  for(i in 1:length(y.at)){
    plot.xy(
      xy.coords(
        list(
          x= helper.x,
          y= rep(y.at[i], 2)
        )),
      type = "l",
      lty = "dashed",
      lwd = 1,
      col= "grey75"
    )
  }
  
  # axis
  # X
  axis(side = 1,
       at = x.at,
       labels = x.lab,
       cex.axis = 1.2,
       tcl=-.3, 
       col = NA, # supress axis line
       col.ticks = axis.col, 
       col.axis =  axis.col, 
       lwd = 2,
       family="sans",
       font = 2,
       col.lab=text.col
  )
  
  title(ylab = "Significant genes", 
        col.lab=text.col, 
        cex.lab=1.2,
        family="sans",
        font = 2)
  
  # Y
  axis(side = 2,
       at = y.at,  
       las = 1,
       labels = y.lab,
       tcl=-.3, 
       col = NA, # supress axis line
       col.ticks = axis.col, 
       col.axis =  axis.col, 
       lwd = 2,
       family="sans",
       font = 2,
       col.lab=text.col
  )

  # draw lines from each origin metabolites to their quotient
  
  for(i in plot.dat[quotient.def == "Quotient", unique(group)]){
    
    # i <- "Q6"
    x1 <- base.plot.dat[metabolite == (i), group.num.jitter]
    y1 <- base.plot.dat[metabolite == (i), sig.genes.log]
    
    for(j in plot.dat[quotient.def != "Quotient" & group == (i), metabolite]){
      
      # readline(prompt="Press [enter] to continue")
      # message(paste0("Drawing connection of metabolite ", j, " to quotient ", i))

      # j <- "Met"
      x2 <- base.plot.dat[metabolite == (j), group.num.jitter]
      y2 <- base.plot.dat[metabolite == (j), sig.genes.log]
      
      # lwd for quot.sig T mehr, F weniger
      # mto = more than one 
      mto <- plot.dat[group==(i) & metabolite == (j), quot.sigs]
      
      # color, lwd and alpha
      col <- plot.dat[group==(i), unique(group.col)]
      
      line.lwd <- ifelse(col!="less", 2.5, ifelse(mto==TRUE, 1.5, 1))

      line.col <- ifelse(col != "less", alpha("indianred4", .6), 
                         ifelse(mto==TRUE, alpha("indianred", .6),
                                alpha("grey45", .5)))
      
      # draw each line
      plot.xy(
        xy.coords(
          list(
            x= c(x1, x2),
            y= c(y1, y2)
          )),
        type = "c",
        lwd = line.lwd,
        col= line.col
      )
    }
  }
  
  
  
  # Invariant Size
  
  # color based on new perc
  
  # pch based on type
  base.plot.dat[, plot.pch := as.numeric(
    as.character(
      factor(
        metab.class, 
        levels = c("Amino Acid", "Acylcarnitine", "Mix"), 
        labels=c(21,22,23)
        )
      )
    )]
  
  # fill/contour color for aa/ac/mix
  base.plot.dat[, dot.col := ifelse(metab.class.q == "aa", "darkgoldenrod",
                                     ifelse(metab.class.q == "ac", "cadetblue4", 
                                            "darkolivegreen4"))]
  
  # get a gradient for the new.perc variable for all quotients but not for the AA/AC only dots - these become grey
  base.plot.dat[quotient.def!="Quotient", quot.fill := "grey95"]
  new.perc.pal <- colorRampPalette(colors = c("snow1", "darkolivegreen4"))
  
  #This adds a column of color values
  n.colors <- base.plot.dat[quotient.def=="Quotient",uniqueN(new.perc)]
  all.cols <- new.perc.pal(n.colors)
  plot.cols <- all.cols[as.numeric(cut(base.plot.dat[quotient.def=="Quotient",new.perc],breaks = n.colors))]
  
  # add the correct colors to the data
  base.plot.dat[quotient.def=="Quotient", quot.fill := plot.cols]
  
  # plot points of each metabolite
  plot.xy(
    xy.coords(
      list(
        x= base.plot.dat$group.num.jitter,
        y= base.plot.dat$sig.genes.log
      )),
    type = "p",
    pch = base.plot.dat$plot.pch, 
    lwd = 2.5,
    col= alpha(base.plot.dat$dot.col, 0.9), 
    bg = alpha(base.plot.dat$quot.fill, 0.9),
    cex = 2
  )
  
  #' ## Add densities to the side!
  
  # plot(x = c(4, 8.5), y = c(0, 1.5), type = "n", 
  #      xlab = "Sepal Length", ylab = "Density", 
  #      main = "Density plot of multiple groups"); 
  # 
  # d <- lapply(split(iris$Sepal.Length, iris$Species), density)
  # Map(function(dens, col) polygon(dens, col = col), 
  #     dens = d, col = c('pink', 'lightgreen', 'lightblue')); 
  # 
  # legend("topright", levels(iris$Species), fill = c('pink', 'lightgreen', 'lightblue'))
  # 
  
  #' ## Add legends
  
  # legend for percentage of unique genes in quotient 
  legend(
    x = "bottomleft",
    title = "In Quotients:",
    title.adj =  0,
    title.col = text.col,
    y.intersp = .9,
    legend = c("0% new genes", 
               "0-50% new genes", 
               "50-100% new genes"),
    pt.cex = 2,
    col = "darkolivegreen4",
    pt.bg = new.perc.pal(3),
    pch = 21,
    pt.lwd = 2.5, 
    cex = 1.2, 
    bty = "n", 
    text.col = text.col, 
    inset = c(0, -.25),
    horiz = F
  )
  
  # legend for metab class color
  legend(
    y.intersp = .9,
    x = "bottomright", 
    legend = c("Amino Acid", 
               "Acylcarnitine", 
               "Mix Quotient", 
               "Amino Acid Quotient", 
               "Acylcarnitine Quotient"), 
    pt.bg = c("Amino Acid" = alpha("grey97", 0.9), 
              "Acylcarnitine" = alpha("grey97", 0.9),
              "Mix Quotient" = alpha("grey97", 0.9),
              "Amino Acid Quotient" = alpha("grey97", 0.9), 
              "Acylcarnitine Quotient" = alpha("grey97", 0.9)), 
    col = c("Amino Acid" = alpha("darkgoldenrod", .8), 
            "Acylcarnitine" = alpha("cadetblue4", .8),
            "Mix Quotient" = alpha("darkolivegreen4", .8),
            "Amino Acid Quotient" = alpha("darkolivegreen4", .8), 
            "Acylcarnitine Quotient" = alpha("darkolivegreen4", .8)), 
    pt.lwd = 2.5, 
    pch = c(21, 22, 23, 21, 22), 
    bty = "n", 
    pt.cex = 2, 
    cex = 1.2, 
    text.col = text.col,
    inset = c(0, -.27),
    horiz = F
    )
  
  # for line type
  legend(
    x = "bottom", 
    title = "In Quotients:",
    title.adj =  0,
    title.col = text.col,
    legend = c("Less Significant Genes", 
               "More Significante Genes than this part",
               "More Significant Genes than all parts"),
    lwd = c(1,
            1.5,
            3),
    col = c(alpha("grey45", .5),
            alpha("indianred", .7), 
            "indianred4"),
    pt.lwd = 1.5, 
    lty = "solid",
    bty = "n", 
    pt.cex = 2, 
    cex = 1.2, 
    text.col = text.col, 
    inset=c(-0.25),
    horiz = F
  )
  
}
