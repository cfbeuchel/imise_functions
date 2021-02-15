rescale_metabolites <- function(unscaled_data, scaled_data, cols){
  
  res <- purrr::map2(.x = unscaled_data[ , .SD, .SDcols = cols],
                     .y = scaled_data[ , .SD, .SDcols = cols],
                     .f = function(x, y) {
                       
                       # multiply the transformed values, now SD = 1, with the SD of the untransformed values to scale them back to its original SD
                       # original sd
                       my.sd <- sd(x, na.rm = T)
                       my.mean <- mean(x, na.rm = T)
                       
                       # transformed col * original sd
                       my.y <- y * my.sd + my.mean
                       return(my.y)
                     })
  
  # is a list, needs to be a dt
  setDT(res)
  return(res)
}