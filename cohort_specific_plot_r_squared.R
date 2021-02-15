plot_r_squared <- function(myData = all.res, confounders, correction, rSquaredCutoff){
  
  # rearrange data for plotting
  molten.r.squared <- melt(data = myData,
                           measure.vars = names(myData)[names(myData) %in% confounders],
                           id.vars = c("cohort", "response"))
  
  # build plot 
  p1 <- ggplot(na.omit(molten.r.squared), aes(x = variable,
                                              y = as.numeric(value),
                                              fill = variable)) +
    # geom_point(size = 0, alpha = 0) + # zeichnet punkte für label und entfernt die wiede
    geom_hline(yintercept = rSquaredCutoff, col = "red") +
    geom_boxplot() +
    ggtitle(label = "Covariate partial r-squared, all cohorts",
            subtitle = paste0("Variable selection based on ", correction, " cutoff")) + 
    # facet_grid(~cohort) +
    # scale_y_log10(breaks = scales::pretty_breaks(10)) +
    theme_carl() +
    theme(axis.text.x = element_text(angle= 45, hjust = 1, size = 11),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  # call plot
  return(p1)
  
}