glm_mediation_loop <- function(y,
                               x1,
                               x2,
                               c,
                               id,
                               testPlan,
                               type = c("binomial", "gaussian")){
  
  stopifnot(type %in% c("binomial", "gaussian"))
  stopifnot(uniqueN(testPlan$pheno)==1)
  
  # get the correct overlap 
  v <- venn4(y[[id]], x2[[id]], x1[[id]], c[[id]],plotte = FALSE)$q1
  
  # match the correct ids
  m1 <- match(v, y[[id]])
  m2 <- match(v, x1[[id]])
  m2.5 <- match(v, x2[[id]])
  m3 <- match(v, c[[id]])
  
  # create single data object
  dat.mod <- cbind(
    y[m1, ],
    x1[m2, ],
    x2[m2.5, ],
    c[m3, ]
  )
  
  
  testPlan[,index:=1:.N]
  res <- lapply(testPlan$index, function(i){
    
    pred.m <- testPlan[index==(i), metabolite]
    pred.g <- testPlan[index==(i), gx.probe]
    
    
    # build the correct model fomula
    formula.covars <- paste(names(c)[names(c) != id] , collapse = " + ")
    formula.response <- paste0("as.numeric( ", names(y)[names(y) != id], " ) ~ ")
    
    # edit 200714
    # old:
    # formula.predictor <- paste(pred.m, " + ", pred.g, " + ")
    formula.predictor <- paste0(pred.m, " + ", pred.g, " + ", pred.g, ":", pred.m, " + ")
    
    mod.formula <- as.formula(paste0(formula.response, formula.predictor, formula.covars))
    
    # get all the IDs of non-na data for the correct LR test by tanking all NAs over all involved variables
    non.na.ids <- !dat.mod[, apply((is.na(.SD)), 1, any), .SDcols=c(id, pred.m, pred.g, names(y)[names(y)!=id], names(c)[names(c)!=id])]
    
    # fit the model and extract coefficients
    mod.glm <- glm(formula = mod.formula, family = type, data = dat.mod)
    mod.coef <- as.data.frame(coefficients(summary(mod.glm)))[, ]
    
    # edit 200714: add interaction effect and extract coefficient
    # old: mod.coef <- mod.coef[rownames(mod.coef) %in% c(pred.m,pred.g),]
    mod.coef <- mod.coef[rownames(mod.coef) %in% c(pred.m,
                                                   pred.g, 
                                                   paste0(pred.m, ":", pred.g),
                                                   paste0(pred.g, ":", pred.m)
                                                   ),]
    
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
      length(mod.glm)
      N <- summary(mod.glm)$df.null
      
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
    
    id1 <- testPlan[index==(i), paste0(metabolite, "__", gx.probe)]
    id2 <- i
    mod.coef[x == pred.g, x.type := "gx.probe"]
    mod.coef[x == pred.m, x.type := "metabolite"]
    mod.coef[x == paste0(pred.m,":",pred.g), x.type := "interaction"]
    
    # annotate the mediating factor
    mod.coef[x == pred.g, m := pred.m]
    mod.coef[x == pred.m, m := pred.g]
    
    mod.coef[, `:=`(
      id.name = id1,
      id.number = id2
    )]
    
  })
  
  # consolidate results
  res <- rbindlist(res)
  res[, y := names(y)[names(y) != id]]
  res <- res[,.SD,.SDcols=c(
    "id.number",
    "id.name",
    "y",
    "x",
    "m",
    "x.type",
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
