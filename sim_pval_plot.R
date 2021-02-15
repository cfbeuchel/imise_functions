sim_pval_plot <- function(data){
  
  # change all.scenarios model label
  new.names <- c("Ranks", 
                 "Asinh-transformation",
                 "Inverse-normal-transformation",
                 "Dichotomization",
                 "Categorization")
  
  names(data)[match(c("spearman.cor.pval.zero.inflated",
                      "linmod.pval.zero.inflated",
                      "linmod.pval.zero.inflated.transformed",
                      "dicho.pval",
                      "logit.pval"),
                    names(data))] <- new.names
  
  # show comparison with full data <- linmod, dicho, spearman
  asm <- melt(data, id.vars = c("effekt",
                                "n",
                                "proz0",
                                "categart",
                                "categnum",
                                "batch.effect",
                                "combat.run",
                                "combat.sub.run",
                                "num"),
              measure.vars = new.names,
              # measure.vars = names(data)[grep(x = names(data),
              #                                 pattern = "pval")],
              variable.name = "model",
              value.name = "pval", F)
  
  # plotting colors
  pc <- c("#46024E", "#460E54", "#461859", "#45205E", "#452763", "#432E68",
          "#41346D", "#3F3A72", "#3B4177", "#37477B", "#314D7F", "#2A5383",
          "#205987", "#105F8A", "#00658D", "#006B90", "#007193", "#007795",
          "#007C97", "#008298", "#008899", "#008D9A", "#00929A", "#00989A",
          "#009D9A", "#00A299", "#00A698", "#00AB96", "#00B095", "#00B492",
          "#00B890", "#00BC8D", "#00C089", "#0AC486", "#31C881", "#47CB7D",
          "#59CF78", "#69D273", "#77D46E", "#85D768", "#93DA62", "#A0DC5C",
          "#ADDE56", "#B9E050", "#C5E14A", "#D1E245", "#DDE340", "#E8E43C",
          "#F3E438", "#FDE333")
  
  # pc <- c("#6D6D6D", "#F9F9F9")
  
  # creat data
  for.plot <- asm[categart == "quantile" &
                    categnum == 10 &
                    batch.effect == 1 &
                    model %in% new.names, ]
  
  
  # remove unwanted effect sizes etc
  # for.plot <- for.plot[
  #   # effekt %in% c(0,0.02, 0.05, 0.1, 0.3) &
  #   model %nin% "Linear Model w/ asinh-transformation"
  #   ]
  
  # power of assoc
  # for.plot[, power := sum(pval <= 0.05)/.N, by = c("effekt", "proz0", "model")]
  model.power <- for.plot[, .(power = round(sum(pval <= 0.05)/.N, 3),
                              xlim = min(pval)), by = c("effekt", "proz0", "model")]
  model.power[ , xlim.f := sort(xlim)[6], by = c("effekt")]
  
  # plot
  pval.density <- ggplot(for.plot, aes(y = model,
                                       x = pval,
                                       fill = ifelse(..x.. <= log10(0.05), "sig", "no sig"))) +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",
                                  # calc_ecdf = T,
                                  # quantiles = c(0.05), 
                                  panel_scaling = T) +
    scale_fill_manual(name = element_blank(),
                      values = c(pc[1], pc[45]),
                      labels = c("p > 0.05", "p < 0.05")) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    # facet_grid(proz0 ~ effekt, scales = "free_x") + 
    facet_grid(proz0 ~ effekt,
               scales = "free_x",
               labeller = label_bquote(rows = .(proz0) ~ "% zeros",
                                       cols = beta ==.(effekt))) +
    xlab(label = "P-value") +
    geom_text(data = model.power,
              aes(x=xlim.f, y = model, label = power),
              size = 3, nudge_y = .4, hjust = 0) +
    theme_carl() +
    theme(
      axis.ticks.x = element_line(),
      axis.title.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "bottom",
      panel.border = element_blank())
}
