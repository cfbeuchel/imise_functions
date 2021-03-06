multivariate_assoc_metab_covar <- function(dometab, docovar, data, r_on_server = F) {
  
  #=================================================================#
  # SETUP ----------------------------------------------------------#
  
  # number of cores used in parallel processing
  if (r_on_server == T) {
    number.of.cores <- 10
  } else if (r_on_server == F) {
    number.of.cores <- parallel::detectCores()/2
  }
  
  # Testplan
  todo <- expand.grid(covar = docovar, metab = dometab, stringsAsFactors = F)
  setDT(todo)
  
  all_metab_covar_multi <- lapply(X = unique(todo$metab), FUN = function(mymetab){
    
    # mymetab <- unique(todo$metab)[8] # debug
    # define metabolite and confounders
    mycovars <- todo[metab == mymetab, covar]
    
    # (generalized) linear Model with tryCatch
    res <- tryCatch(
      {
          #================================================================#
          # LINEAR MODEL --------------------------------------------------#
          
          # define the formula 
          myformula.lm <- as.formula(paste0(mymetab," ~ ", paste(mycovars, collapse = "+")))
          
          # Fit the linear model 
          mod.raw <- lm(myformula.lm, data = data)
          mod.vif <- as.data.table(vif(mod.raw), keep.rownames = T)
          mod <- summary(mod.raw)
          
          # Clean the output
          # get the estimate from the raw model output because it includes all missings!
          coeffs.raw.lm <- as.data.frame(mod.raw$coefficients[grep(pattern = paste0(mycovars, collapse = "|"), x = names(mod.raw$coefficients), value = F)])
          # coeffs.raw.lm <- as.data.frame(mod.raw$coefficients[names(mod.raw$coefficients) %in% mycovars])
          setDT(coeffs.raw.lm, keep.rownames = T)
          setnames(coeffs.raw.lm, c("rn", "estimate"))
          
          # get the other statistiscs from the tidy output
          # coeffs.lm <- mod$coefficients[rownames(mod$coefficients) %in% mycovars, ]
          coeffs.lm <- mod$coefficients[grep(pattern = paste0(mycovars, collapse = "|"), x = names(mod.raw$coefficients), value = F), ]
          covar.order <- rownames(mod$coefficients)[grep(pattern = paste0(mycovars, collapse = "|"), x = names(mod.raw$coefficients), value = F)]
          
          # stupid special case with only one covar 
          if(length(mycovars) == 1) {
            coeffs.raw.lm[ , `:=` (std.error.lm = coeffs.lm[which(names(coeffs.lm) == "Std. Error")],
                                   pval.lm = coeffs.lm[which(names(coeffs.lm) == "Pr(>|t|)")])]
          } else {
            
            coeffs.lm <- as.data.frame(coeffs.lm)
            setDT(coeffs.lm, keep.rownames = T)
            # match the coeffs coeffs.raw.lm
            coeffs.raw.lm[ , `:=` (std.error.lm = coeffs.lm[match(coeffs.raw.lm$rn, covar.order) , `Std. Error`],
                                   pval.lm = coeffs.lm[match(coeffs.raw.lm$rn, covar.order) , `Pr(>|t|)`])]
          }
          
          # remove double mentions from coeffs.raw.lm
          # check which coeffs.raw.lm$rn is mentioned more than once in coeffs.raw.lm and take the one with the smalles p-value
          exclude.rows <- sapply(mycovars, function(i){
            db <- grep(pattern = i, x = coeffs.raw.lm$rn)
            dbi <- coeffs.raw.lm[db, which(pval.lm == min(pval.lm))]
            exclude <- db[db %nin% db[dbi]]
          })
          
          # get all row indices to exclude
          exclude.rows <- unlist(exclude.rows)
          
          # filter rosw
          coeffs.raw.lm <- coeffs.raw.lm[!(exclude.rows), ]
          
          # change rownames of vif dt in case something is off again
          if(ncol(mod.vif) == 2){
            names(mod.vif) <- c("rn", "GVIF")
          }
          
          # Enter results into dt #  HIER BAUSTELLE!!!!
          res <- data.table::data.table(metab = mymetab,
                                        term = mycovars,
                                        estimate = coeffs.raw.lm$estimate,
                                        std.error = coeffs.raw.lm$std.error.lm,
                                        r.squared = mod$adj.r.squared,
                                        p.value = coeffs.raw.lm$pval.lm,
                                        vif = mod.vif[match(rn, mycovars), GVIF],
                                        n = length(mod$residuals),
                                        comment = NA)
          
          # END OF: LINEAR MODEL ------------------------------------------#
          #================================================================#
      },
      error = function(cond){
        res <- data.table(metab = mymetab,
                          term = mycovars,
                          estimate = NA,
                          std.error = NA,
                          r.squared = NA,
                          p.value = NA, 
                          vif = NA,
                          n = NA,
                          comment = stringr::str_replace_all(as.character(cond), "\n", "  "))
        return(res)
      }
    ) # End of tryCatch()
    return(res)
  })
  
  all_metab_covar_multi <- rbindlist(all_metab_covar_multi)
  return(all_metab_covar_multi)
}
