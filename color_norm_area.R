color_dens_area <- function(density_function, area_below, ...) {
  
  # get coordinates on x, simple line
  cord.x <- c(-3, seq(-3, area_below, 0.01), area_below)
  
  # coordinates on y need to follow the specified density function
  cord.y <- c(0, eval(density_function)(seq(-3, area_below, 0.01), ...), 0)
  
  # draw the curve
  curve(eval(density_function)(x, ...), xlim = c(-3, 3))
  
  # draw the colored polygon marking the area below the input value
  polygon(cord.x, cord.y, col = 'skyblue')
}



color_dens_area(density_function = pgamma, area_below =  2, shape = 1.2, rate = 0.5)
