compare_limma_matrixeqtl <- function(
  datLimma,
  datMatrix,
  checkCohort,
  checkMetabolite,
  checkStatistic = c("pvalue", "beta")
  ){
  
  # change
  if(checkStatistic == "pvalue") {
    
    # in case of pval compare, rename the cols
    setnames(datMatrix, "pvalue", "toCompare", skip_absent = TRUE) 
    setnames(datLimma, "P.Value", "toCompare", skip_absent = TRUE) 
    my.plot <- qqplot
    
    
  } else {
    
    # in case of beta, rename beta cols
    setnames(datMatrix, "beta", "toCompare", skip_absent = TRUE) 
    setnames(datLimma, "logFC", "toCompare", skip_absent = TRUE) 
    my.plot <- plot
    
  }
  
  allDuplicatedEntries(datLimma[cohort == (checkCohort) & var == (checkMetabolite), ilmn])
  order.limma <- match(datMatrix[cohort == (checkCohort) & snps == (checkMetabolite), gene],
                       datLimma[cohort == (checkCohort) & var == (checkMetabolite), ilmn])
  datMatrix[cohort == (checkCohort) & snps == (checkMetabolite), ]
  datLimma[cohort == (checkCohort) & var == (checkMetabolite), ][order.limma, ]
  
  # check
  all(datMatrix[cohort == (checkCohort) & snps == (checkMetabolite), gene] ==
        datLimma[cohort == (checkCohort) & var == (checkMetabolite), ][order.limma, ilmn]
  )
  
  # check some metabolites
  my.plot(x = datMatrix[cohort == (checkCohort) & snps == (checkMetabolite), toCompare],
         y = datLimma[cohort == (checkCohort) & var == (checkMetabolite), ][order.limma, toCompare], 
         xlab = paste0(checkCohort, " metabolite ", checkMetabolite, " (MatrixEQTL)"), 
         ylab = paste0(checkCohort, " metabolite ", checkMetabolite, " (LIMMA)"), 
         main = paste0(checkCohort, ": ", checkMetabolite, " ", checkStatistic))

  # change
  if(checkStatistic == "pvalue") {
    
    # in case of pval compare, rename the cols
    setnames(datMatrix, "toCompare", "pvalue", skip_absent = TRUE) 
    setnames(datLimma, "toCompare", "P.Value", skip_absent = TRUE) 
    
  } else {
    
    # in case of beta, rename beta cols
    setnames(datMatrix, "toCompare", "beta", skip_absent = TRUE) 
    setnames(datLimma, "toCompare", "logFC", skip_absent = TRUE) 
    
  }
}
