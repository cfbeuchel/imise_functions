sim_pval_plot_selected <- function(data, effect, zeroInflation){
  
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
  
  # colors for plotting -> plot in B/W for journal!
  # plotting colors
  # pc <- c("#262626", "#2A2A2A", "#2F2F2F", "#333333", "#383838", "#3C3C3C", 
  #         "#414141", "#464646", "#4A4A4A", "#4F4F4F", "#545454", "#585858", 
  #         "#5D5D5D", "#626262", "#666666", "#6B6B6B", "#707070", "#747474", 
  #         "#797979", "#7E7E7E", "#828282", "#878787", "#8C8C8C", "#909090", 
  #         "#959595", "#999999", "#9E9E9E", "#A2A2A2", "#A7A7A7", "#ABABAB", 
  #         "#B0B0B0", "#B4B4B4", "#B8B8B8", "#BCBCBC", "#C0C0C0", "#C4C4C4", 
  #         "#C8C8C8", "#CCCCCC", "#D0D0D0", "#D4D4D4", "#D8D8D8", "#DBDBDB", 
  #         "#DEDEDE", "#E2E2E2", "#E5E5E5", "#E8E8E8", "#EBEBEB", "#EDEDED", 
  #         "#EFEFEF", "#F1F1F1")
  # color palette
  pc <- c("#46024E", "#460E54", "#461859", "#45205E", "#452763", "#432E68",
          "#41346D", "#3F3A72", "#3B4177", "#37477B", "#314D7F", "#2A5383",
          "#205987", "#105F8A", "#00658D", "#006B90", "#007193", "#007795",
          "#007C97", "#008298", "#008899", "#008D9A", "#00929A", "#00989A",
          "#009D9A", "#00A299", "#00A698", "#00AB96", "#00B095", "#00B492",
          "#00B890", "#00BC8D", "#00C089", "#0AC486", "#31C881", "#47CB7D",
          "#59CF78", "#69D273", "#77D46E", "#85D768", "#93DA62", "#A0DC5C",
          "#ADDE56", "#B9E050", "#C5E14A", "#D1E245", "#DDE340", "#E8E43C",
          "#F3E438", "#FDE333")
  
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
  
  # select only necessary data
  for.plot.2 <- asm[categart == "quantile" &
                      categnum == 10 & 
                      effekt %in% (effect) &
                      proz0 %in% (zeroInflation) &
                      batch.effect == 1 &
                      model %in% new.names, ]
  
  for.plot.2 <- for.plot.2[
    model %nin% "Linear Model w/\nasinh-transformation"
    ]
  
  # power for fill
  model.power.2 <- for.plot.2[, .(power = round(sum(pval <= 0.05)/.N, 3),
                                  xlim = min(pval)), , by = c("effekt", "proz0", "model")]
  model.power.2[ , xlim.f := sort(xlim)[2], by = c("effekt")]
  
  # reduced plot
  pval.density.2 <- ggplot(for.plot.2,
                           aes(y = model,
                               x = pval, 
                               fill = ifelse(..x.. <= log10(0.05), "sig", "no sig"))) + 
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",
                                  panel_scaling = T) +
    geom_text(data = model.power.2,
              aes(x = xlim.f, y = model, label = power),
              size = 4, nudge_y = .2, hjust = 0) +
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
    theme_carl() +
    theme(text = element_text(size=15),
          axis.ticks.x = element_line(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "bottom",
          panel.border = element_blank())
  
}
