theme_gxMetab <- function(
  bg.col = "grey100",
  axis.col = "grey30",
  text.col = "grey29",
  strip.bg.col = "grey98",
  grid.col = "grey30"
  ){
  
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(), 
    panel.border         = element_blank(),
    panel.background     = element_rect(fill = bg.col,
                                        color=NA),
    panel.grid.major     = element_line(color = grid.col,
                                        linetype = "dotted"),
    axis.ticks           = element_line(color = axis.col, 
                                        size=.6),
    axis.ticks.length    = unit(.20,"cm"),
    axis.text            = element_text(family = "sans",
                                        colour = text.col,
                                        face = "bold"),
    axis.title           = element_text(family = "sans",
                                        colour = text.col,
                                        face = "bold"),
    strip.background     = element_rect(color=NA,
                                        fill=strip.bg.col),
    strip.text           = element_text(family = "sans",
                                        colour = text.col,
                                        face = "bold"),
    legend.justification = "top"
  )
}
