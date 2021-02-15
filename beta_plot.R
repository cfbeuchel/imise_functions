beta_plot <- function(x, y, xQval, yQval, fdrControl = 0.05, showPlot = TRUE){
  
  # checks prior to plot
  stopifnot({
    !anyNA(
      c(x, 
        y, 
        xQval, 
        yQval, 
        fdrControl)
    )
    length(x) == length(y)
    length(xQval) == length(yQval)
    length(x) == length(xQval)
    length(y) == length(yQval)
    min(xQval) >= 0 & max(xQval) <= 1
    min(yQval) >= 0 & max(yQval) <= 1
    fdrControl >= 0 & fdrControl <= 1
  })
  
  # consolidate for plotting and annotate significancy based on fdrControl for color plotting
  dat <- data.table(
    x = x,
    y = y,
    x.q = xQval,
    y.q = yQval,
    coloring = ifelse((xQval <= fdrControl) & (yQval <= fdrControl), "x & y",
                      ifelse(xQval <= fdrControl, "x",
                             ifelse(yQval <= fdrControl, "y", "none")
                      )
    )
  )  
  
  # check for zero significant associations in x & y -> nothing to plot
  if(all(dat$x.q >= fdrControl) & all(dat$y.q >= fdrControl)){
    stop(paste0("No significant associations found for q<=", 
                fdrControl, 
                ". Maybe try setting a higher FDR?"))
  }
  
  # calculate r-squared of model in assocs that are at least significant in x or y
  # mod <- dat[coloring != "none", summary(lm(x~y))]
  mod <- dat[x.q < 0.05 | y.q < 0.05, summary(lm(x~y))]
  cor.pearson <- dat[x.q < 0.05 | y.q < 0.05, cor.test(x,y,method="pearson")]
  
  # get the annotation
  r2 <- bquote(R^2* '='~.(signif(mod$r.squared,2)))
  
  # create plot - uses default theme 
  # use theme_set(theme_light(base_size = 12, base_family = "Helvetica")) or similar
  #p <- ggplot(dat[coloring != "none", ], aes(x = x,
  p <- ggplot(dat[x.q < 0.05 | y.q < 0.05, ], aes(x = x,
                       y = y,
                       col = coloring,
                       alpha=coloring)) + 
    scale_alpha_manual(values=c("x"=0.8,"y"=0.8,"x & y"=0.8,"none"=.1),
                       guide=F) +
    geom_point(size = 2) + 
    geom_hline(yintercept = 0, col  = "grey55") + 
    geom_vline(xintercept = 0, col  = "grey55") + 
    geom_rug(length=unit(0.02,"npc"),
             outside=FALSE,
             sides="tr") +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    geom_abline(intercept = 0, slope = 1, lty = 2, col = "#7C0000") + 
    annotate("text", 
             x = -Inf, 
             y = Inf,
             hjust = -0.3,
             vjust = 2,
             label = deparse(r2),
             size = 6, 
             parse = T) +
    guides(col = guide_legend(
        title = bquote("q" ~ ""<="" ~ .(fdrControl) ~ "in:"))) +
    theme(legend.position="right",
          plot.margin = margin(1, 1, 1, 1, "cm"))

  
  if(showPlot == TRUE){
    print(p)
  }
  
  return(list(
    plot = p, # returns plot for easier modification with labels etc
    data = dat, # returns dat used for producing the plot
    model = mod, # for getting the exact R2 or something else produced by the model
    cor = cor.pearson
  ))
}

