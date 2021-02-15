glm_loop <- function(y,x,c,id, type = c("binomial", "gaussian")){
  
  stopifnot(type %in% c("binomial", "gaussian"))
  
  # get the correct overlap 
  v <- venn3(y[[id]], x[[id]], c[[id]],plotte = FALSE)$q1
  
  # match the correct ids
  m1 <- match(v, y[[id]])
  m2 <- match(v, x[[id]])
  m3 <- match(v, c[[id]])
  
  # create single data object
  dat.mod <- cbind(
    y[m1, ],
    x[m2, ],
    c[m3, ]
  )
  
  res <- lapply(names(x)[names(x) != id], function(i){
    
    # build the correct model fomula
    formula.covars <- paste(names(c)[names(c) != id] , collapse = " + ")
    formula.response <- paste0("as.numeric( ", names(y)[names(y) != id], " ) ~ ")
    formula.predictor <- paste(i, " + ")
    mod.formula <- as.formula(paste0(formula.response, formula.predictor, formula.covars))
    
    # get all the IDs of non-na data for the correct LR test by tanking all NAs over all involved variables
    non.na.ids <- !dat.mod[, apply((is.na(.SD)), 1, any), .SDcols=c(id, i, names(y)[names(y)!=id], names(c)[names(c)!=id])]
    
    # fit the model and extract coefficients
    if(type=="gaussian"){
      link.type <- gaussian()
    } else if(type=="binomial"){
      # link.type <- binomial(link = "log")
      link.type <- binomial(link = "logit")
    }
    
    mod.coef <- tryCatch({
      
      # try to fit the model
      mod.glm <- glm(formula = mod.formula, family = link.type, data = dat.mod)
      mod.coef <- as.data.frame(coefficients(summary(mod.glm)))[i, ]
      setDT(mod.coef, keep.rownames = TRUE)
      setnames(mod.coef, c("x", "beta", "beta.se", "test.statistic", "p.value"))
      mod.coef[, n := length(residuals(mod.glm))]
      
      # get different r2 depending on model type
      if (type=="gaussian") {
        
        # get the correct r2
        mod.sum <- summary(mod.glm)
        r2 <- with(mod.sum, 1 - deviance/null.deviance)
        r2.adj <- with(mod.sum, 1 - (1 - r2) * ((df.null)/(df.residual)))
        mod.coef[,`:=`(
          r.squared.1 = r2,
          r.squared.2 = r2.adj,
          r.squared.3 = beta^2 / ( beta^2 + n * beta.se^2 ),
          r.squared.type = "r2;r2.adj;arnd"
        )]
        
      } else {
        
        # estimate the NULL LogLikelihood
        null.formula <- as.formula(paste0(formula.response, " 1"))
        mod.null <- glm(formula = null.formula, family = type, data = dat.mod[non.na.ids, ])
        
        # get a pseudo r2
        K <- length(c(names(c)[names(c) != id], i))
        # length(mod.glm)
        N <- nrow(mod.glm$data)
        
        # methods("logLik")
        # getAnywhere("logLik.lm")
        L1 <- logLik(mod.glm)
        L0 <- logLik(mod.null)
        LR <- lmtest::lrtest(mod.glm,mod.null)$Chisq[2]
        
        # get the pseudo-r2
        r2.mcfadden <- 1 - ((L1)/L0)
        r2.mcfadden.adj <- (1 - (( L1 - K )/L0))[1]
        r2.nagelkerke <- (1-exp(-LR/N))/(1-exp(-(-2*L0)/N))[1]
        mod.coef[,`:=`(
          r.squared.1 = r2.mcfadden.adj,
          r.squared.2 = r2.nagelkerke,
          r.squared.3 = beta^2 / ( beta^2 + n * beta.se^2 ),
          r.squared.type = "mcfadden.adj;nagelkerke;arnd"
        )]
        
      }
      
      # empty comment for later merging
      mod.coef[, comment := NA]
      
      # return successful glm result
      return(mod.coef)
      
    }, error = function(cond){
      mod.coef <- data.table(
        x = i,
        beta = NA,
        beta.se = NA,
        test.statistic = NA,
        p.value = NA,
        n = NA,
        r.squared.1 = NA,
        r.squared.2 = NA,
        r.squared.3 = NA,
        r.squared.type = NA,
        comment = gsub(x = as.character(cond),pattern =  "\n", 
                       replacement = " ",
                       fixed = TRUE))
      
      # return failed glm result
      return(mod.coef)
    })
    
  })
  
  # consolidate results
  res <- rbindlist(res)
  res[, y := names(y)[names(y) != id]]
  res <- res[,.SD,.SDcols=c("y",
                            "x",
                            "beta",
                            "beta.se",
                            "test.statistic",
                            "p.value",
                            "r.squared.1",
                            "r.squared.2",
                            "r.squared.3",
                            "r.squared.type",
                            "n")]
  
  return(res)
  
}
