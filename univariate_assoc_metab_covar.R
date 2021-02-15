univariate_assoc_metab_covar <- function(data, dometab, docovar, r_on_server = F) {
  
  # number of cores used in parallel processing
  if (r_on_server == T) {
    number.of.cores <- 10
  } else if (r_on_server == F) {
    number.of.cores <- parallel::detectCores()/2
  }
  
  # Testplan
  todo <- expand.grid(covar = docovar, metab = dometab, stringsAsFactors = F)
  setDT(todo)
  
  # First check, whether single covariables are significantly associated with the metabolites
  
  #================================================================#
  # SINGLE COVARIABLES --------------------------------------------#
  all_metab_covar <- parallel::mclapply(unique(todo$metab), function(mymetab){
    
    # mymetab <- unique(todo$metab)[8] # debug
    # create toto-list based on first mymetab
    subtodo <- todo[metab == mymetab, ]
    
    all_covar <- lapply(unique(subtodo$covar), function(mycovar){
      
      # mycovar <- unique(subtodo$covar)[2]
      
      res <- tryCatch(
        {
          #================================================================#
          # LINEAR MODEL --------------------------------------------------#
          
          # Fit the linear model 
          myformula.lm <- as.formula(paste0(mymetab, " ~ ", mycovar))
          mod <- summary(lm(formula = myformula.lm, data = data))
          
          # Clean the output
          coeffs.lm <- as.data.frame(mod$coefficients)
          setDT(coeffs.lm, keep.rownames = T)
          estimate.lm <- unlist(coeffs.lm[rn %in% (grep(pattern = mycovar, coeffs.lm$rn, value = T)), "Estimate"])
          std.error.lm <- unlist(coeffs.lm[rn %in% (grep(pattern = mycovar, coeffs.lm$rn, value = T)), "Std. Error"])
          pval.lm <- unlist(coeffs.lm[rn %in% (grep(pattern = mycovar, coeffs.lm$rn, value = T)), "Pr(>|t|)"])
          
          # Enter results into dt
          res <- data.table::data.table(metab = mymetab,
                                        term = mycovar,
                                        estimate = estimate.lm,
                                        std.error = std.error.lm,
                                        r.squared = mod$r.squared,
                                        p.value = pval.lm,
                                        n = length(mod$residuals),
                                        comment = NA)
          
          # END OF: LINEAR MODEL ------------------------------------------#
          #================================================================#
        },
        error = function(cond){
          res <- data.table(metab = mymetab,
                            term = mycovar,
                            estimate = NA,
                            std.error = NA,
                            r.squared = NA,
                            p.value = NA,
                            n = NA,
                            comment = stringr::str_replace_all(as.character(cond), "\n", "  "))
          return(res)
        }
      ) # End of tryCatch()
      
      # in case of factor, this will be longer than 1 -> only use the most important factor
      if(nrow(res) != 1){

        # get the result with the highest r²
        res <- res[which(p.value == min(p.value, na.rm = T)), ]
      }
      
      # return res to all_covar as a result of the lapply loop
      return(res)
      
    }) # End of covar-lapply
    
    # bind together the output of one metabolite and all covars 
    all_covar2 <- rbindlist(all_covar, use.names = T)
    
    # this is the result that will be collected in all_metab_covar
    return(all_covar2)
    
  }, mc.cores = number.of.cores, mc.cleanup = T) # End of metab-lapply mc: , mc.cores = number.of.cores, mc.cleanup = T
  
  # END OF: SINGLE COVARIABLES ------------------------------------#
  #================================================================#
  
  # Bind the final data.table together with all metabolites
  all_metab_covar <- rbindlist(all_metab_covar, use.names = T)
  
  #output
  return(all_metab_covar)
}
