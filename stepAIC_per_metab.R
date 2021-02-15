stepAIC_per_metab <- function(metab) {
  
  # debug: metab <- "Q35"
  
  # get all covariates for each cohort
  all.covars <- list(a1 = annot$covar.annot[a1.in.confounder.analysis == T, covariate[covariate %in% reshaped.data$term]],
                     b3 = annot$covar.annot[b3.in.confounder.analysis == T, covariate[covariate %in% reshaped.data$term]],
                     sorb = annot$covar.annot[sorb.in.confounder.analysis == T, covariate[covariate %in% reshaped.data$term]])
  
  # create formulas for each model
  all.formulas <- sapply(all.covars, function(x){
    as.formula(paste(metab, " ~ ", paste(x, collapse = " + ")))
  })
  
  # create upper (max model) and lower (min model) part of scope for stepAIC
  upper <- all.formulas
  lower <- list(a1 = as.formula(paste(metab," ~ 1")),
                b3 = as.formula(paste(metab, " ~ 1")),
                sorb = as.formula(paste(metab, " ~ 1")))
  
  # create scope variable
  all.scopes <- map2(.x = upper, 
                     .y = lower,
                     .f = function(x,y){
                       
                       # create scope for all three cohorts
                       list(upper = x,
                            lower = y)
                     })
  
  # set cohort identifier
  all.cohorts <- list(a1 = "a1",
                      b3 = "b3",
                      sorb = "sorb")
  
  # fit the models
  # fit the upper models
  all.mods <- pmap(.l = list(my.form = all.formulas,
                             my.covar = all.covars,
                             my.cohort = all.cohorts),
                   
                   .f = function(my.form,
                                 my.covar,
                                 my.cohort,
                                 my.data = metab.covar,
                                 my.metab = metab){
                     
                     # debug
                     # my.form <- all.formulas$a1
                     # my.covar <- all.covars$a1
                     # my.cohort <- all.cohorts$a1
                     # my.data <- metab.covar
                     
                     # col definition needed for ommiting NAs
                     my.cols <- c(my.metab, my.covar)
                     
                     # model fit without NAs
                     my.mod <- lm(formula = my.form, data = my.data[cohort == my.cohort, na.omit(.SD), .SDcols = my.cols])
                     return(my.mod)
                   })
  
  # run stepaic
  all.best.mods <- pmap(.l = list(my.cohort = all.cohorts,
                                  my.mod = all.mods,
                                  my.scope = all.scopes,
                                  my.covar = all.covars),
                        .f = function(my.cohort,
                                      my.mod,
                                      my.scope,
                                      my.covar,
                                      my.data = metab.covar,
                                      my.metab = metab) {
                          
                          # debug
                          # my.cohort <- all.cohorts$a1
                          # my.mod <- all.mods$a1
                          # my.covar <- all.covars$a1
                          # my.scope <- all.scopes$a1
                          # my.data <- metab.covar
                          
                          # reassign my.cols because stepAIC is confused by .SDcols = <something>
                          my.cols <- c(my.metab, my.covar)
                          
                          # run forward and backwards selection stepAIC
                          invisible(my.best.mod <- (MASS::stepAIC(object = my.mod, scope = my.scope, direction = "backward")))
                          my.res <- list(metab = my.metab,
                                         covars = attr(my.best.mod$terms, "term.labels"),
                                         AIC = AIC(my.best.mod))
                          return(my.res)
                        }
  )
  return(all.best.mods)
}