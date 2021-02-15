var_na_heatmaps <- function(dat) {
  
  na.res <- dat[ , lapply(.SD, function(x){
    sum(is.na(x))/length(x)}), by = "FileName", .SDcols = batch.metabolites]
  
  # plot Nas over all metabolites and batches
  na.res <- na.res[, .SD, .SDcols = 2:length(names(na.res))]
  my.na.matrix <- na.res[, as.matrix(.SD, dimnames = list(unique(dat$FileName),names(na.res)))]
  rownames(my.na.matrix) <- unique(dat$FileName)
  heatmap(1 - my.na.matrix, distfun = NA, Rowv = NA, Colv = NA, scale = "none", main = "NAs in % per batch")
  
  # check the variance of my data, since 0 variance causes ComBat to not converge
  # calculate the variance of each metabolite per batch
  var.res <- dat[ , lapply(.SD, function(x){
    var(x)}), by = FileName,  .SDcols = batch.metabolites]
  
  # heatmap of variance 
  var.res <- var.res[, .SD, .SDcols = 2:length(names(var.res))]
  my.var.matrix <- var.res[, as.matrix(.SD, dimnames = list(unique(dat$FileName),names(na.res)))]
  rownames(my.var.matrix) <- unique(dat$FileName)
  heatmap(my.var.matrix, distfun = NA, Rowv = NA, Colv = NA, scale = "row", main = "Metabolite variance by batch")
}