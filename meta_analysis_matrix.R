meta_analysis <- function(
  phenotype,
  assocData
){
  
  ### Compute server --------------------------------------------------------------------------------
  computer <- "forostar" # specify exactly, i. e. with 'MRO' suffix (if applicable)
  
  ### Path for QC data ------------------------------------------------------------------------------
  pathProjectQC <- "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1512_meta_meta/04_ge_gx_metab/190808_gx_metab_Association/"
  

  ### Methods to be performed -----------------------------------------------------------------------
  Fisher <- TRUE
  ZScore <- TRUE
  FEM    <- TRUE
  REM    <- TRUE
  
  genomicControl <- FALSE # only, if lambda_gc > 1
  
  ### Number format in final data table -------------------------------------------------------------
  signDig <- 6 # number of significant digits
  
  ### Phenotype name as part of meta data file name to be saved
  # phenotype <- "..." # is assigned in function input
  
  ### Additional suffix for meta data file name to be saved
  tag <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  
  ###################################################################################################
  ### Data import -----------------------------------------------------------------------------------
  
  # split into list
  assocData <- assocData[metabolite %in% (phenotype), ]
  data <- lapply(unique(assocData$cohort), function(i){
    
    obj <- assocData[cohort == (i), .(markerID, codedAll, noncodedAll, eaf, imputed, infoscore, n, beta, se, p)]
    return(obj)
  })
  
  # Unite "markerID"
  allMarkerID <- c()
  for (i in 1 : length(data)) {
    
    allMarkerID <- c(allMarkerID, as.character(data[[i]][, markerID]))
  }
  
  # Match data by "allMarkerID"
  allMarkerID <- as.factor(unique(allMarkerID))
  for (i in 1 : length(data)) {
    
    match     <- match(allMarkerID, data[[i]][, markerID])
    data[[i]] <- data[[i]][match]
  }
  allMarkerID <- as.character(allMarkerID)
  
  # >>>>>>>>>> IMPORTANT! <<<<<<<<<<
  nameStudies <- unique(assocData$cohort)

  k <- length(data) # total number of studies
  m <- nrow(data[[1]]) # number of matched SNPs
  ###################################################################################################
  
  
  
  
  
  ###################################################################################################
  ### Collect single data ---------------------------------------------------------------------------
  # Generate data frame with markerID, total sample size, n-weighted MAF, minimum MAF, imputation indicator, n-weighted info score, and minimum info score
  n         <- c()
  maf       <- c()
  imputed   <- c()
  infoscore <- c()
  for (i in 1 : k) {
    
    n                 <- cbind(n, data[[i]][, n])
    
    dum               <- data[[i]][, eaf]
    dum1              <- dum > 0.5
    dum1[is.na(dum1)] <- FALSE
    dum[dum1]         <- 1 - dum[dum1]
    maf               <- cbind(maf, dum)
    
    imputed   <- cbind(imputed, data[[i]][, imputed])
    infoscore <- cbind(infoscore, data[[i]][, infoscore])
  }
  numberStudies      <- apply(!is.na(n), 1, sum)
  totalN             <- apply(n, 1, sum, na.rm = TRUE)
  nWeight            <- n / totalN
  nWeightedMAF       <- apply((nWeight * maf), 1, sum, na.rm = TRUE)
  minMAF             <- apply(maf, 1, min, na.rm = TRUE)
  nWeightedInfoScore <- apply((nWeight * infoscore), 1, sum, na.rm = TRUE)
  minInfoScore       <- apply(infoscore, 1, min, na.rm = TRUE)
  
  metaData <- data.table(markerID = allMarkerID,
                         numberStudies = numberStudies,
                         totalN = totalN,
                         nWeightedMAF = nWeightedMAF,
                         minMAF = minMAF,
                         nWeightedInfoScore = nWeightedInfoScore,
                         minInfoScore = minInfoScore)
  
  # Study results
  for (i in 1 : k) {
    
    metaData <- data.table(metaData, signif(data[[i]][, beta], digits = signDig), signif(data[[i]][, se], digits = signDig), signif(data[[i]][, p], digits = signDig), data[[i]][, n], signif(data[[i]][, eaf], digits = signDig), signif(maf[, i], digits = signDig), imputed[, i], signif(infoscore[, i], digits = signDig))
  }
  
  # Study specific column names
  colNamesStudyResults <- c()
  for (i in 1 : k) {
    
    colNamesStudyResults <- c(colNamesStudyResults, 
                              paste0("beta.", nameStudies[i]), 
                              paste0("se.", nameStudies[i]), 
                              paste0("p.", nameStudies[i]), 
                              paste0("n.", nameStudies[i]), 
                              paste0("eaf.", nameStudies[i]), 
                              paste0("maf.", nameStudies[i]), 
                              paste0("imputed.", nameStudies[i]), 
                              paste0("infoscore.", nameStudies[i]))
  }
  setnames(metaData, c((ncol(metaData) - 8 * k + 1) : ncol(metaData)), colNamesStudyResults)
  ###################################################################################################
  
  ###################################################################################################
  ### Fisher's method -------------------------------------------------------------------------------
  if (Fisher == TRUE) {
    
    p <- c()
    for (i in 1 : k) {
      
      p <- cbind(p, (-2 * log(data[[i]][, p])))
    }
    chisq <- apply(p, 1, sum, na.rm = TRUE)
    
    pFisher <- pchisq(chisq, 2 * k, lower.tail = FALSE)
  }
  
  
  
  ### Z score ---------------------------------------------------------------------------------------
  if (ZScore == TRUE) {
    
    num <- c()
    den <- c()
    for (i in 1 : k) {
      
      num <- cbind(num, qnorm(1 - data[[i]][, p] / 2) * sign(data[[i]][, beta]) * sqrt(data[[i]][, n]))
      den <- cbind(den, data[[i]][, n])
    }
    num <- apply(num, 1, sum, na.rm = TRUE)
    den <- apply(den, 1, sum, na.rm = TRUE)
    den <- sqrt(den)
    Z   <- num / den
    
    pZ <- 2 * pnorm(abs(Z), lower.tail = FALSE)
  }
  
  
  
  ### Fixed effects model ---------------------------------------------------------------------------
  if (FEM == TRUE) {
    
    theta <- c()
    w     <- c()
    
    for (i in 1 : k) {
      
      theta0 <- data[[i]][, beta]
      se0    <- data[[i]][, se]
      w0     <- 1 / se0^2 # das ist IVW
      
      if (genomicControl == TRUE) {
        
        lambdaGC <- median((theta0 / se0)^2, na.rm = TRUE) / 0.456
        if (lambdaGC > 1) {
          
          w0 <- (1 / (lambdaGC * se0^2))
        }
      }
      
      theta <- cbind(theta, theta0)
      w     <- cbind(w, w0)
    }
    
    thetaFEM      <- apply((w * theta), 1, sum, na.rm = TRUE) / apply(w, 1, sum, na.rm = TRUE)
    # To avoid numerical problems for calculation of 'w * theta / w = theta', 'thetaFEM = theta' is forced
    # -->
    dum           <- numberStudies == 1
    thetaFEM[dum] <- apply(theta[dum,],1,  sum, na.rm = T)
    # <--
    seFEM         <- sqrt(1 / apply(w, 1, sum, na.rm = TRUE))
    
    ZFEM <- thetaFEM / seFEM
    
    pFEM <- 2 * pnorm(abs(ZFEM), lower.tail = FALSE)
    
    # Additionally, Cochran's Q statistic and I^2 statistic will be calculated for data output
    Q       <- apply((w * (theta - thetaFEM)^2), 1, sum, na.rm = TRUE)
    I2      <- (Q - (numberStudies - 1)) / Q
    dum     <- (I2 < 0) | is.na(I2)
    I2[dum] <- 0
  }
  
  
  
  ### Random effects model --------------------------------------------------------------------------
  if (REM == TRUE) {
    
    # START: from FEM (without Genomic Control)
    theta <- c()
    w     <- c()
    for (i in 1 : k) {
      
      theta <- cbind(theta, data[[i]][, beta])
      w     <- cbind(w, (1 / data[[i]][, se]^2))
    }
    thetaFEM <- apply((w * theta), 1, sum, na.rm = TRUE) / apply(w, 1, sum, na.rm = TRUE)
    # To avoid numerical problems for calculation of 'w * theta / w = theta', 'thetaFEM = theta' is forced
    # -->
    dum           <- numberStudies == 1
    thetaFEM[dum] <- apply(theta[dum,],1,  sum, na.rm = T)
    # <--
    # END: from FEM (without Genomic Control)
    
    Q         <- apply((w * (theta - thetaFEM)^2), 1, sum, na.rm = TRUE)
    tau2      <- (Q - (numberStudies - 1)) / (apply(w, 1, sum, na.rm = TRUE) - (apply((w^2), 1, sum, na.rm = TRUE) / apply(w, 1, sum, na.rm = TRUE)))
    dum       <- (Q < (numberStudies - 1)) | is.na(tau2)
    tau2[dum] <- 0
    I2        <- (Q - (numberStudies - 1)) / Q
    dum       <- (I2 < 0) | is.na(I2)
    I2[dum]   <- 0
    
    wREM <- matrix(nrow = m, ncol = k)
    for (j in 1 : k) {
      
      for (i in 1 : m) {
        
        wREM[i, j] <- 1 / (1 / w[i, j] + tau2[i])
      }
    }
    
    thetaREM <- apply((wREM * theta), 1, sum, na.rm = TRUE) / apply(wREM, 1, sum, na.rm = TRUE)
    seREM    <- sqrt(1 / apply(wREM, 1, sum, na.rm = TRUE))
    
    ZREM <- thetaREM / seREM
    
    pREM <- 2 * pnorm(abs(ZREM), lower.tail = FALSE)
  }
  ###################################################################################################
  
  
  
  
  
  ###################################################################################################
  ### Collect meta data -----------------------------------------------------------------------------
  # Meta results
  if (Fisher == TRUE) { # add p value
    
    metaData$pFisher <- signif(pFisher, digits = signDig)
  }
  
  if (ZScore == TRUE) { # add Z score and p value
    
    metaData$Z  <- signif(Z, digits = signDig)
    metaData$pZ <- signif(pZ, digits = signDig)
  }
  
  if (FEM == TRUE) { # add beta, standard error, p value, Cochran's Q statistic and its p value, and I^2 statistic
    
    metaData$betaFEM    <- signif(thetaFEM, digits = signDig)
    metaData$seFEM      <- signif(seFEM, digits = signDig)
    metaData$pFEM       <- signif(pFEM, digits = signDig)
    metaData$CochransQ  <- signif(Q, digits = signDig)
    metaData$pCochransQ <- signif(pchisq(Q, df = (numberStudies - 1), lower.tail = FALSE), digits = signDig)
    metaData$I2         <- signif(I2, digits = signDig)
  }
  
  if (REM == TRUE) { # add beta, standard error, p value, ...
    
    metaData$betaREM <- signif(thetaREM, digits = signDig)
    metaData$seREM   <- signif(seREM, digits = signDig)
    metaData$pREM    <- signif(pREM, digits = signDig)
    
    if (FEM == FALSE) { # ...and - if not already done - Cochran's Q statistic, its p value, and I^2 statistic
      
      metaData$CochransQ  <- signif(Q, digits = signDig)
      metaData$pCochransQ <- signif(pchisq(Q, df = (numberStudies - 1), lower.tail = FALSE), digits = signDig)
      metaData$I2         <- signif(I2, digits = signDig)
    }
  }
  
  # Save as R data object and file in 'Results' folder
  # assign(phenotype, metaData)
  metaData[, metabolite := (phenotype)]
  return(metaData)
  
  
  # pathProjectResults <- paste0(pathProjectQC, "/../Results/") # path for results folder
  # dir.create(pathProjectResults) # create results folder
  # fn                 <- paste0(pathProjectResults, "GWASMA_", phenotype, "_", tag) # file name for meta data to be saved
  # # with 'phenotype' and 'tag'
  # 
  # save(list = phenotype, file = paste0(fn, ".RData"))
  # fwrite(metaData, file = paste0(fn, ".txt"), na = "NA")
  # gzip(paste0(fn, ".txt"), destname = paste0(fn, ".gz"))
  ###################################################################################################
  
}
