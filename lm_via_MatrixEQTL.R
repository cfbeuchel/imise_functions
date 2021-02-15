lm_via_MatrixEQTL =  function(x,
                              y,
                              conf= NULL,
                              idcolumn = "id",
                              max_sclice_x = Inf,
                              max_sclice_y = Inf,
                              max_sclice_conf = Inf,
                              max_pval=1,
                              showplot=T,
                              stoponcheck = T) {
  # idcolumn = "id"; max_sclice_x = Inf; max_sclice_y = Inf; max_sclice_conf = Inf;max_pval=1
  

  # minimal example
  # source("R:/genstat/02_projekte/1703_ifb_lifea1_lifeb3_sorbs_sugar/lm_via_MatrixEQTL.R")
  # 
  # ge = data.table(id = paste0("ilmn_",letters[1:5]), matrix(rnorm(50), ncol = 10  ))
  # ge
  # 
  # metab = data.table(id = paste0("metabo_",letters[1:4]), matrix(rnorm(40), ncol = 10  ))
  # metab
  # 
  # covars = data.table(id = paste0("covs_",letters[1:3]), matrix(rnorm(30), ncol = 10  ))
  # covars
  # 
  # # mit adjustierung 
  # assoc = lm_via_MatrixEQTL(x = ge, y = metab, conf = covars, idcolumn = "id", stoponcheck = T)
  # assoc
  # 
  # 
  # # ohne adjustierung
  # assoc = lm_via_MatrixEQTL(x = ge, y = metab,idcolumn = "id", stoponcheck = T)
  # assoc
  
  makeMetrischSlicedDaten = function(x, idcolumnsub = idcolumn) {
    library(data.table)
    setDT(x)
    library(MatrixEQTL)
    metrische_var_indinfo = as.matrix(x[,-idcolumnsub, with = F])
    row.names(metrische_var_indinfo) = x[,get(idcolumnsub)]
    hh(metrische_var_indinfo)
    sliced_snps =  SlicedData$new()
    sliced_snps = sliced_snps$CreateFromMatrix( metrische_var_indinfo)
    sliced_snps
  }
  
  
  library(data.table)
  
  stopifnot(ncol(x)==ncol(y))
  if(is.null(conf) ==F) stopifnot(ncol(x)==ncol(conf))
  
  
  setDT(x)
  setDT(y)
  setDT(conf)
  
  setcolorder(x, c(idcolumn, sort(setdiff(names(x), idcolumn))))
  setcolorder(y, c(idcolumn, sort(setdiff(names(y), idcolumn))))
  stopifnot(identical(names(y), names(x)))
  
  if(nrow(conf) >0)  {
    setcolorder(conf, c(idcolumn, sort(setdiff(names(conf), idcolumn))))
    stopifnot(identical(names(conf), names(x)))
    }
  
  
  sliced_x = makeMetrischSlicedDaten(x)
  sliced_y = makeMetrischSlicedDaten(y)
  if(nrow(conf) >0) sliced_conf = makeMetrischSlicedDaten(conf) else sliced_conf = SlicedData$new()
  
  if(ncol(x) > max_sclice_x) {
    sliced_x$ResliceCombined(sliceSize = max_sclice_x)  
  }
  
  if(ncol(y) > max_sclice_y) {
    sliced_y$ResliceCombined(sliceSize = max_sclice_y)  
  }
  
  if(is.null(conf) ==F) {
    if(ncol(conf) > max_sclice_conf) {
      sliced_conf$ResliceCombined(sliceSize = max_sclice_conf)  
    }
  } 
  
  
  
  assoc = Matrix_eQTL_engine(
    snps = sliced_x,
    gene = sliced_y,
    cvrt = sliced_conf, 
    output_file_name = NULL, 
    pvOutputThreshold = max_pval,
    useModel = modelLINEAR,
    noFDRsaveMemory=F,
    pvalue.hist = "qqplot")
  
  if(showplot) plot(assoc)
  assoc_pvals = assoc$all$eqtls
  setDT(assoc_pvals)
  
  assoc_pvals[,beta_se := beta / statistic]
  assoc_pvals[,dfFull := assoc$param$dfFull]
  assoc_pvals[,r2 := (statistic / sqrt( dfFull + statistic^2 ))^2]
  
  assoc_pvals
  
  
  ## checks
  if(nrow(conf) >0) {
    df_vgl1 = data.table(x = as.numeric(t(x[1,-idcolumn, with = F])),
                         y = as.numeric(t(y[1,-idcolumn, with = F])),
                         t(conf[,-idcolumn, with = F])
    )} else  {df_vgl1 = data.table(x = as.numeric(t(x[1,-idcolumn, with = F])),
                                   y = as.numeric(t(y[1,-idcolumn, with = F])),
                                   conf = 1)}
  df_vgl1
  
  
  vgl1_zeile = paste0('summary(lm( y ~  x +', paste(setdiff(names(df_vgl1), c("x", "y")), collapse= " + "), ' , data=df_vgl1))')
  vgl1 = eval(parse(text = vgl1_zeile ))
  vgl1 = vgl1$coefficients
  vgl1
  x_oriname1 = as.character(x[1,idcolumn, with = F])
  y_oriname1 = as.character(y[1,idcolumn, with = F])
  
  check1 = identical(round(vgl1["x", "t value"],5), assoc_pvals[snps ==x_oriname1 & gene == y_oriname1, round(statistic,5)])
  statcol = grep("value", colnames(vgl1), value = T)
 
  
  if(nrow(conf) >0) {
    df_vgl2 = data.table(x = as.numeric(t(x[nrow(x),-idcolumn, with = F])),
                         y = as.numeric(t(y[nrow(y),-idcolumn, with = F])),
                         t(conf[,-idcolumn, with = F])
    )} else  {df_vgl2 = data.table(x = as.numeric(t(x[nrow(x),-idcolumn, with = F])),
                                   y = as.numeric(t(y[nrow(y),-idcolumn, with = F])),
                                   conf = 1)}
  df_vgl2
  
  
  vgl2_zeile = paste0('summary(lm( y ~  x +', paste(setdiff(names(df_vgl2), c("x", "y")), collapse= " + "), ' , data=df_vgl2))')
  vgl2 = eval(parse(text = vgl2_zeile ))
  vgl2 = vgl2$coefficients
  vgl2
  x_oriname2 = as.character(x[nrow(x),idcolumn, with = F])
  y_oriname2 = as.character(y[nrow(y),idcolumn, with = F])
  
  check2 = identical(round(vgl2["x", "t value"],5), assoc_pvals[snps ==x_oriname2 & gene == y_oriname2, round(statistic,5)])
  statcol = grep("value", colnames(vgl2), value = T)
 
  
  if(check1 & check2) message("For two examples (",x_oriname1, " vs. ", y_oriname2, " and ", x_oriname2, " vs.", y_oriname1, ") R ANOVA and MatrixEQTL ANOVA incl. covariate are identical. Great.") else {
    message("For two examples (",x_oriname1, " vs. ", y_oriname2, " and ", x_oriname2, " vs.", y_oriname1, ") R ANOVA and MatrixEQTL ANOVA incl. covariate were NOT identical. BAD.")
    print(vgl1)
    print(assoc_pvals[snps ==x_oriname1 & gene == y_oriname1, ])
    print("and\n")
    print(vgl2)
    print(assoc_pvals[snps ==x_oriname2 & gene == y_oriname2, ])
    
    
    }
    if(stoponcheck) stopifnot(check1 & (is.na(vgl1[,statcol][1])==F))
    if(stoponcheck) stopifnot(check2 & (is.na(vgl2[,statcol][1])==F))
    
  return(assoc_pvals)
  
} 
