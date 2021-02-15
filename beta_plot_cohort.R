beta_plot <- function(x, y, xQval, yQval, cohort, fdrControl = 0.05, showPlot = TRUE){
  
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
    ),
    cohort = cohort
  )  
  
  # check for zero significant associations in x & y -> nothing to plot
  if(all(dat$x.q >= fdrControl) & all(dat$y.q >= fdrControl)){
    stop(paste0("No significant associations found for q<=", 
                fdrControl, 
                ". Maybe try setting a higher FDR?"))
  }
  
  # calculate r-squared of model in assocs that are at least significant in x or y
  mod <- dat[coloring != "none", summary(lm(x~y))]
  
  # get the annotation
  r2 <- bquote(R^2* '='~.(signif(mod$r.squared,2)))
  
  # create plot - uses default theme 
  # use theme_set(theme_light(base_size = 12, base_family = "Helvetica")) or similar
  p <- ggplot(dat[coloring != "none", ], aes(x = x,
                                             y = y,
                                             col = coloring)) + 
    facet_grid(~cohort) + 
    geom_point(alpha = 0.4, size = 2) + 
    geom_hline(yintercept = 0, col  = "grey55") + 
    geom_vline(xintercept = 0, col  = "grey55") + 
    geom_abline(intercept = 0, slope = 1, lty = 2, col = "#7C0000") + 
    annotate("text", 
             x = -Inf, 
             y = Inf,
             hjust = -0.3,
             vjust = 2,
             label = deparse(r2),
             size = 6, 
             parse = T) +
    guides(col = guide_legend(title = bquote("q" ~ ""<="" ~ .(fdrControl) ~ "in:")))
  
  if(showPlot == TRUE){
    print(p)
  }
  
  return(list(
    plot = p, # returns plot for easier modification with labels etc
    data = dat, # returns dat used for producing the plot
    model = mod # for getting the exact R2 or something else produced by the model
  ))
}

