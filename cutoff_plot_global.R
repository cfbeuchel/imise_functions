cutoff_plot_global <- function(plotData){
  
# build plot
cutoff.plot.global <- ggplot(plotData,
                             aes(x = r.squared.cutoff, y = covars)) + 
  geom_point(aes(col = step, pch = step), size = 3, alpha = 0.8) +
  # geom_vline(aes(xintercept = r.squared.cutoff), lwd = .2, col = "red", lty = "dotted") + 
  # geom_hline(aes(yintercept = 10), lwd = .5, col = "black", lty = "dashed") + 
  geom_line(aes(group = step, col = step)) + 
  facet_grid(~correction, labeller = labeller(correction = labels.correction)) + 
  scale_colour_manual(values = plotting.pal.g, labels = c("Univariable", "Multivariable")) +
  scale_shape(labels = c("Univariable", "Multivariable")) +
  scale_y_continuous(breaks = pretty_breaks(10)) + 
  # guides(shape =T) + 
  theme_carl() +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle= 0, hjust = 1),
        axis.line.y = element_blank(),
        legend.title = element_blank(), 
        panel.border = element_blank(),
        # axis.title.x = element_blank(),
        legend.position = "bottom") +
  xlab(label = "R² cutoff") + 
  ylab(label = "Covariates matching filter criteria")

return(cutoff.plot.global)
}
