zero_inflation_percentage <- function(metabolites){
  unlist(metabolites[, lapply(.SD, function(x) sum(x == 0, na.rm = T)/length(x))])
}
