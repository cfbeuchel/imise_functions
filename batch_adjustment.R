batch_adjustment <- function(dat, batch, parameters, batch.date = NA, phenos = NA, my.par.prior = F) {
  
  # designate transposed matrix for adjustment
  my.edata <- dat[, t(as.matrix(.SD)), .SDcols = parameters]
  
  # set batch as character string
  my.batch <- dat[, unname(unlist(.SD)), .SDcols = batch]
  
  # create formula in case phenos are not NA
  if(!all(is.na(phenos))) {
    
    # set phenotypes to adjust
    my.pheno <- dat[, .SD, .SDcols = c(batch, phenos)]
    
    # add covars
    my.formula <- as.formula(paste0("~ 1 + ", paste(phenos, collapse = " + ")))
    
  } else {
    
    # set phenotypes to adjust
    my.pheno <- dat[, .SD, .SDcols = batch]
    
    # only intercept term, no covariates
    my.formula <- ~1
  }
  
  # create model matrix
  my.modcombat <- model.matrix(my.formula, data = my.pheno)
  
  if(!is.na(batch.date)) {
    
    #=================================================================#
    # ANOVA - DATE #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # run an ANOVA with the date
    # batch.date <- substr(dat[, unname(unlist(.SD)), .SDcols = batch], 1, 8)
    
    # loop through all parameters
    all.aov.results <- sapply(parameters, function(my.param) {
      
      # custom formula for each parameter
      my.formula.2 <- as.formula(paste0(my.param, " ~ ", batch.date))
      
      # calculate anova
      my.aov <- anova(aov(my.formula.2, data = dat))
      
      # clean output
      clean.aov <- as.data.table(broom::tidy(my.aov))
      aov.results <- clean.aov[term == "date", p.value]
      
      # return p.vals
      return(aov.results)
    }) # end ANOVA sapply
    
    # summary ANOVA results
    aov.res <- data.table(parameter = names(all.aov.results), p.values = all.aov.results)
    
    # END: ANOVA - DATE #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #=================================================================#
    
  } else {
    aov.res <- NULL
  } # End ANOVA

    if(!all(is.na(phenos))) {
    
    # batch.date <- substr(dat[, unname(unlist(.SD)), .SDcols = batch], 1, 8)
    
    # loop through all parameters
    all.aov.results <- lapply(parameters, function(my.param) {
      
      # custom formula for each parameter
      my.formula.2 <- as.formula(paste0(my.param, " ~ ", paste(phenos, collapse=" + ")))
      
      # calculate anova
      my.aov <- summary(lm(my.formula.2, data = dat))
      
      # clean output
      clean.aov <- as.data.table(broom::tidy(my.aov))
      aov.results <- clean.aov[-1, ]
      aov.results$metab <- my.param
      
      # return p.vals
      return(aov.results)
    }) # end ANOVA sapply
    
    # summary ANOVA results
    aov.pheno.res <- rbindlist(all.aov.results)

  } else {
    aov.pheno.res <- NULL
  } # End ANOVA
  
  #=================================================================#
  # ANOVA - BATCH, BEFORE #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # check the batch effect before adjustment 
  
  # loop through all parameters
  before.aov.results <- sapply(parameters, function(my.param) {
    
    # custom formula for each parameter
    my.formula.3 <- as.formula(paste0(my.param, " ~ ", batch))
    
    # calculate anova
    my.aov <- anova(aov(my.formula.3, data = dat))
    
    # clean output
    clean.aov <- as.data.table(broom::tidy(my.aov))
    single.before.aov.results <- clean.aov[term == (batch), p.value]
    
    # return p.vals
    return(single.before.aov.results)
  }) # end ANOVA sapply
  
  # summary ANOVA results
  before.aov.res <- data.table(parameter = names(before.aov.results), p.values = before.aov.results)
  
  # END: ANOVA - BATCH, BEFORE #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #=================================================================#
  
  # run ComBat for parametric plot
  try(ComBat(dat = my.edata, batch = my.batch, mod = my.modcombat, par.prior = T, prior.plots = TRUE))
  
  # reset plot pars
  par(mfrow = c(1,1))
  
  # run combat non-parametrically
  combat.edata <- ComBat(dat = my.edata, batch = my.batch, mod = my.modcombat, par.prior = my.par.prior, prior.plots = TRUE)
  
  # reset plot pars
  par(mfrow = c(1,1))
  
  # re-enter adjusted data into dat
  dat.dummy <- copy(dat)
  dat.dummy[, (parameters) := as.data.table(t(combat.edata))]
  
  #=================================================================#
  # ANOVA - BATCH, AFTER #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # check the batch effect after adjustment 
  
  # loop through all parameters
  after.aov.results <- sapply(parameters, function(my.param) {
    
    # custom formula for each parameter
    my.formula.3 <- as.formula(paste0(my.param, " ~ ", batch))
    
    # calculate anova
    my.aov <- anova(aov(my.formula.3, data = dat.dummy))
    
    # clean output
    clean.aov <- as.data.table(broom::tidy(my.aov))
    single.after.aov.results <- clean.aov[term == (batch), p.value]
    
    # return p.vals
    return(single.after.aov.results)
  }) # end ANOVA sapply
  
  # summary ANOVA results
  after.aov.res <- data.table(parameter = names(after.aov.results), p.values = after.aov.results)
  
  # END: ANOVA - BATCH, AFTER #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #=================================================================#
  
  if(!all(is.na(phenos))) {
    
    # batch.date <- substr(dat[, unname(unlist(.SD)), .SDcols = batch], 1, 8)
    
    # loop through all parameters
    all.aov.pheno.after <- lapply(parameters, function(my.param) {
      
      # custom formula for each parameter
      my.formula.2 <- as.formula(paste0(my.param, " ~ ", paste(phenos, collapse=" + ")))
      
      # calculate anova
      my.aov <- summary(lm(my.formula.2, data = dat.dummy))
      
      # clean output
      clean.aov <- as.data.table(broom::tidy(my.aov))
      aov.results <- clean.aov[-1, ]
      aov.results$metab <- my.param
      
      # return p.vals
      return(aov.results)
    }) # end ANOVA sapply
    
    # summary ANOVA results
    aov.pheno.res.after <- rbindlist(all.aov.pheno.after)
    
  } else {
    aov.pheno.res.after <- NULL
  } # End ANOVA
  
  
    return(list("ANOVA_results_w_date" = aov.res,
                "batch_adjusted_data" = dat.dummy,
                "ANOVA_before_ComBat" = before.aov.res,
                "ANOVA_after_ComBat" = after.aov.res,
                "LM_covariates_before" = aov.pheno.res,
                "LM_covariates_after" = aov.pheno.res.after
                ))
}
