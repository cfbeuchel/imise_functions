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
    r.squared.one.metab <- sapply(all.covars, function(my.covar){
      
      # my.metab <- all.metabs[1] # debug 
      # my.covar <- all.covars[1] # debug 
      
      # start by building the full and minus one model
      my.full.formula <- as.formula(paste0(my.metab, " ~ ", paste(all.covars, collapse = " + ")))
      my.reduced.formula <- as.formula(paste0(my.metab, " ~ 1",
                                              ifelse(length(all.covars) == 1, "", " + "),
                                              paste(all.covars[which(all.covars != my.covar)], collapse = " + ")))
      
      # fit the two models
      my.full.model <- summary(lm(formula = my.full.formula, data = data[cohort == my.cohort, ]))
      my.reduced.model <- summary(lm(formula = my.reduced.formula, data = data[cohort == my.cohort, ]))
      
      # extract the explained variance 
      my.full.r.squared <- my.full.model$r.squared
      my.full.r.squared.adj <- my.full.model$adj.r.squared
      my.reduced.r.squared <- my.reduced.model$r.squared
      
      # substract r-squares to get explained variance of my.covar
      my.r.squared <- my.full.r.squared - my.reduced.r.squared
      
      res <- list(covar.r2 = my.r.squared,
                  model.r2 = my.full.r.squared.adj)
      
      # return value
      return(res)
    }, USE.NAMES = T) # end of inner loop over all predictors
    
    # sperate and joint together
    covar.r2 <- unlist(r.squared.one.metab["covar.r2",])
    model.r2 <- unlist(unique(r.squared.one.metab["model.r2",]))
    
    # join together
    covar.r2 <- c(covar.r2, total.model = model.r2)
    
    # format for output, one col per predictor, one row per response
    my.output <- data.table(t(covar.r2))
    # my.output <- data.table(t(r.squared.one.metab))
    
    if(verbose == T) {
      # return to parent function
      message(paste(my.metab), " is done!")
    }
    
    return(my.output)
  }) # end of middle loop over all responses
  
  # format output
  my.output <- rbindlist(r.squared.all.metabs)
  
  # add a metabolite column
  my.output[ , `:=`(response = all.metabs,
                    cohort = my.cohort,
                    r.squared.cutoff = rSquaredCutoff)]
  
  # return cohort output
  return(my.output)
  
} # end of function
