cutoff_plot <- function(plotData){
  
  # build plot
  cutoff.plot <- ggplot(plotData, aes(x = cutoff, y = covars)) + 
    geom_point(aes(fill = factor(cohort), pch = factor(cohort), size = factor(cohort)), alpha = 0.8) +
    # geom_vline(aes(xintercept = cutoff), lwd = .2, col = "red", lty = "dotted") + 
    geom_line(aes(group = factor(cohort), col = factor(cohort), lty = factor(cohort))) + 
    # facet_wrap(step~correction, labeller = labeller(step = labels.assoc, correction = labels.correction)) + 
    facet_grid(~correction, labeller = labeller(correction = labels.correction)) +
    theme_carl() +
    scale_y_continuous(breaks = pretty_breaks(10)) +
    scale_colour_manual(values = rep("black", 4), labels = c("LIFE-Adult", "LIFE-Heart", "Sorb cohort", "Combined")) +
    scale_fill_manual(values = rev(plotting.pal), labels = c("LIFE-Adult", "LIFE-Heart", "Sorb cohort", "Combined")) +
    scale_shape_manual(values = c(21,22,23,24), labels = c("LIFE-Adult", "LIFE-Heart", "Sorb cohort", "Combined")) +
    # scale_shape_manual(values = c(15,16,17,18), labels = c("LIFE-Adult", "LIFE-Heart", "Sorb cohort", "Combined")) +
    scale_size_manual(values = c(3,3,3,5), labels = c("LIFE-Adult", "LIFE-Heart", "Sorb cohort", "Combined")) + 
    scale_linetype_manual(values = c("longdash", "longdash", "longdash", "solid"), labels = c("LIFE-Adult", "LIFE-Heart", "Sorb cohort", "Combined")) +
    # guides(size = F) +
    theme_carl() +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle = 0, hjust = 0),
          legend.title = element_blank(), # axis.title.x = element_blank(),
          legend.position = "bottom",
          panel.border = element_blank(),
          axis.line.y = element_blank()) +
    xlab(label = "(partial-)r² cutoff") + 
    ylab(label = "Covariates matching filter criteria")
  
  # return plot
  return(cutoff.plot)
}

