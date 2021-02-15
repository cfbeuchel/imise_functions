inverse_normal_transform <- function(data, cols){
  
  # Loop -------------------------------
  transformed.cols <- data[, lapply(.SD, function(x){
    
    # x <- combat.data$y.0infl_6  # debug
    one.transformed.col <- bestNormalize::orderNorm(x)$x.t
    
    # return to "transformed.cols"-output
    return(one.transformed.col)
    
    # end of lapply-loop
  }), .SDcols = cols]
  # End: Loop --------------------------
  
  # do not change source dt
  data.dummy <- copy(data)
    
  # enter new transformed cols
  data.dummy[ , names(transformed.cols) := transformed.cols]
  
  return(data.dummy)
}