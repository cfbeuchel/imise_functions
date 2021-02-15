correlate_meta_statistics <- function(assocDat, 
                                      statistic,
                                      # controlFDR=0.05,
                                      useCores=4,
                                      verbose=TRUE){
  
  # get all metabolite combinations
  test.plan <- 
    expand.grid(m1=unique(assocDat$metabolite),
                m2=unique(assocDat$metabolite),
                stringsAsFactors = FALSE)
  setDT(test.plan)
  test.plan[,index:=1:.N]
  
  # debug, combination with only 2 observations! cor.test will fail
  # test.plan <- test.plan[4999,] # 
  
  index.seq <- round(seq(1,nrow(test.plan),length.out=100),0)
  index.max <- nrow(test.plan)
  
  # effect size correlations? by gene or metabolite
  all.comb <- parallel::mclapply(test.plan$index, function(x){
    
    # get the metabolites
    i1 <- test.plan[index==(x),m1]
    i2 <- test.plan[index==(x),m2]
    
    # data
    dat.mod <- data.table(
      x = assocDat[numberStudies>1 & metabolite == (i1), unlist(.SD), .SDcols=statistic],
      y = assocDat[numberStudies>1 & metabolite == (i2), unlist(.SD), .SDcols=statistic],
      x.q = assocDat[numberStudies>1 & metabolite == (i1), unlist(.SD), .SDcols="significant"],
      y.q = assocDat[numberStudies>1 & metabolite == (i2), unlist(.SD), .SDcols="significant"]
    )
    
    if(x %in% index.seq & verbose==TRUE){
      print(paste0("Done: ", which(x==index.seq), "%"))
    }
    
    one.combination <- tryCatch({
      
      # Use lm instead and extract beta and r2
      mod <- dat.mod[x.q==TRUE|y.q==TRUE, summary(lm(x~y))]
      cor.mod <- dat.mod[x.q==TRUE|y.q==TRUE, cor.test(x,y,method="pearson")]
      
      # get beta, r2 and pvalue
      # enter into data
      one.combination <- data.table(
        m1=(i1),
        m2=(i2),
        rho=cor.mod$estimate,
        p.rho=cor.mod$p.value,
        beta = mod$coefficients[2,1],
        r2 = mod$r.squared,
        signed.r2 = mod$r.squared * sign(mod$coefficients[2,1]),
        n=length(mod$residuals),
        p.lm = mod$coefficients[2,4],
        comment=NA
      )
      
      return(one.combination)
    },
    error=function(cond){
      
      # create empty object for erroneous objects
      one.combination <- data.table(
        m1=i1,
        m2=i2,
        rho=NA,
        p.rho=NA,
        beta=NA,
        r2=NA,
        signed.r2=NA,
        n=NA,
        p.lm=NA,
        comment=gsub(x = as.character(cond),pattern =  "\n",replacement = "  "))
      
      return(one.combination)
      
    }) # end tryCatch
    
    # return to lapply
    return(one.combination)

    # 
    # if(
    #   dat.mod[x.q<=(controlFDR)|y.q<=(controlFDR), .N] == 0
    # ){
    #   
    #   # empty entry
    #   one.combination <- data.table(
    #     m1=i1,
    #     m2=i2,
    #     rho=NA,
    #     p.rho=NA,
    #     beta=NA,
    #     r2=NA,
    #     signed.r2=NA,
    #     n=NA,
    #     p.lm=NA
    #   )
    #   
    #   # return empty object
    #   return(one.combination)
    #   
    # }else{
    #   
    #   # Use lm instead and extract beta and r2
    #   mod <- dat.mod[x.q<=(controlFDR)|y.q<=(controlFDR), summary(lm(x~y))]
    #   cor.mod <- dat.mod[x.q<=(controlFDR)|y.q<=(controlFDR), cor.test(x,y,method="pearson")]
    #   
    #   # get beta, r2 and pvalue
    #   # enter into data
    #   one.combination <- data.table(
    #     m1=(i1),
    #     m2=(i2),
    #     rho=cor.mod$estimate,
    #     p.rho=cor.mod$p.value,
    #     beta = mod$coefficients[2,1],
    #     r2 = mod$r.squared,
    #     signed.r2 = mod$r.squared * sign(mod$coefficients[2,1]),
    #     n=length(mod$residuals),
    #     p.lm = mod$coefficients[2,4]
    #   )
    #   return(one.combination)
    # } # end else
    
  },mc.cores = useCores, mc.cleanup = T) # end lapply
  
  all.combinations <- rbindlist(all.comb)
  all.combinations[,q.rho:=p.adjust(p.rho, method = "BH")]
  all.combinations[,q.lm:=p.adjust(p.lm, method = "BH")]
  
  return(all.combinations)
  
}
