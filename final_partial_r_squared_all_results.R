final_partial_r_squared <- function(cohort,
                                    responses,
                                    predictors,
                                    data,
                                    rSquaredCutoff,
                                    verbose = T) {
  
  # debug/example
  # cohort <- "a1"
  # predictors <- reshaped.data[!is.na(a1.bonferroni), term]
  # responses <- unique(my.data$metab)
  # data <- metab.covar
  
  # iterativer vorgang: erst r2 bilden, dann das klein
  
  
  # rename for historical compatibility
  my.cohort <- cohort
  all.covars <- predictors
  all.metabs <- responses
  
  # middle loop over all responses, LHS etc
  r.squared.all.metabs <- lapply(all.metabs, function(my.metab) {
    
    # inner loop over all predictors
    r.squared.one.metab <- lapply(all.covars, function(my.covar){
      
      # my.metab <- all.metabs[1] # debug 
      # my.covar <- all.covars[1] # debug 
      
      # start by building the full and minus one model
      my.full.formula <- as.formula(paste0(my.metab, " ~ ", paste(all.covars, collapse = " + ")))
      my.reduced.formula <- as.formula(
        paste0(
          my.metab, " ~ 1",
          ifelse(length(all.covars) == 1, "", " + "),
          paste(all.covars[which(all.covars != my.covar)], collapse = " + ")))
      
      # fit the two models
      my.full.model <- summary(lm(formula = my.full.formula, data = data[cohort == my.cohort, ]))
      my.reduced.model <- summary(lm(formula = my.reduced.formula, data = data[cohort == my.cohort, ]))
      
      # extract the explained variance 
      my.full.r.squared <- my.full.model$r.squared
      my.full.r.squared.adj <- my.full.model$adj.r.squared
      my.reduced.r.squared <- my.reduced.model$r.squared
      
      # clean output
      res <- broom::tidy(my.full.model)
      setDT(res)
      
      # in case of categorical, select the most significant one
      res <- res[term %in% grep(pattern = (paste0("^",my.covar, "[0-9]$")), x = res$term, value = T), ]
      res <- res[p.value == min(p.value),]
      
      # substract r-squares to get explained variance of my.covar
      my.r.squared <- my.full.r.squared - my.reduced.r.squared
      
      # consolidate output
      res[, `:=`(cohort = my.cohort,
                 metab = my.metab,
                 term = my.covar,
                 term.r.squared = my.r.squared,
                 model.r.squared = my.full.r.squared.adj,
                 n = length(my.full.model$residuals))]
      
      # return value
      return(res)
    }) # end of inner loop over all predictors
    
    my.output <- rbindlist(r.squared.one.metab)
    
    if(verbose == T) {
      # return to parent function
      message(paste(my.metab), " is done!")
    }
    
    return(my.output)
  }) # end of middle loop over all responses
  
  # format output
  my.output <- rbindlist(r.squared.all.metabs)
  
  # column order for easier readability
  setcolorder(x = my.output, 
              neworder = c("cohort", "metab", "term", "estimate", "std.error", "statistic", "p.value",  
                            "term.r.squared", "model.r.squared", "n"))
  
  # add a metabolite column
  my.output[ , `:=`(r.squared.cutoff = rSquaredCutoff)]
  
  # return cohort output
  return(my.output)
  
} # end of function
