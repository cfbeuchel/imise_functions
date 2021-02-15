sorb_metabolites_relationship_adjust <- function(dat, metabolites) {
  
  ###
  # Start preparing the data for the relationship adjustment
  ###
  
  # checke, das kein wert kleiner -10000, da das als NA wert verwendet wird, added 21.12.18
  mincheck <- min(unlist(dat[,..metabolites]), na.rm = T)
  if(mincheck <= -10000) stop("Adapt NA codes - this is -10000 but minimum value of data is ", mincheck)
  
  
  # Load population info
  pop <- read.table("/net/ifs1/san_projekte/projekte/genstat/02_projekte/1107_sorben_boettcher/population.txt", header = TRUE)
  pop_id <- as.character(pop$SampleID[pop$Sample == 'Sorben'])
  pop_id <- as.numeric(pop_id)
  
  # load relatedness matrix
  vw <- read.table("/net/ifs1/san_projekte/projekte/genstat/02_projekte/1107_sorben_boettcher/sorben_verwandtschaft.txt", header = FALSE)
  
  # ???
  N <- length(vw[1,])
  vw <- vw + diag(1, N)
  vw <- vw + t(vw) - diag(1, N)
  e <- eigen(vw)
  d <- e$values
  # tail(d)
  d[e$values < 0] <- min(e$values[e$values > 0])
  eig <- e$vectors %*% diag(d) %*% t(e$vectors)
  # eig[1:5, 1:5]
  vw <- eig
  
  # Load identifier gene data #
  ind <- read.table("/net/ifs1/san_projekte/projekte/genstat/04_archiv/sorben/sorben_na27/pheno.txt", header = TRUE)
  d <- !(ind$Populatio == "Sorben500kLeipzig" | ind$Population == "Sorben500kBerlin")
  ind1 <- ind[d, ] #1000k
  d <- !(ind$Population == "Sorben1000k")
  ind2 <- ind[d, ] #500k
  table(ind$SampleID == pop_id) #ok
  ind$SampleID <- as.character(ind$SampleID)
  
  ###
  # prepare the metabolites
  ###
  
  adjust_me <- dat
  # dim(adjust_me)
  
  # Rearrange the data
  filter <- is.element(adjust_me$pid, as.character(ind$SampleID))
  adjust_me <- adjust_me[filter, ]
  # dim(adjust_me)
  
  #  edit Carl 20181223: restrict to ids in data to adjust
  # vw also needs to be restricted
  vw.filter <- pop_id %in% adjust_me$pid
  vw <- vw[vw.filter, vw.filter]
  pop_id <- pop_id[pop_id %in% adjust_me$pid]
  
  # Adjust the order to match the relationship matrix
  Np <- length(pop_id) 
  dum <- data.frame(barcode = adjust_me$pid, Np = 1:Np)
  permu <- c()
  for (i in 1:Np){
    permu[i] <- dum[dum[, 1] == pop_id[i], 2]
  }
  adjust_me <- adjust_me[permu, ]
  
  # Checks
  table(adjust_me$pid == pop_id)       # ok
  
  #  edit Carl 20181223: restrict to ids in data to adjust
  ind <- ind[ind$SampleID %in% pop_id, ]
  
  table(adjust_me$pid == ind$SampleID) # ok
  adjust_me[, pop := ind$Population]
  
  # Ethnic outliers #####################
  ethnic_outlier <- c(2908592, 4556464, 4556481)
  
  # Settings ############################
  zz <- adjust_me
  
  # Select only the metabolites in your data.table
  metab <- zz[, ..metabolites]
  metab[, names(metab) := lapply(.SD, as.double)]
  
  # set NAs to -1s
  metab[, names(metab) := lapply(.SD, function(x){
    x[is.na(x)] <- -10000
    x
  })]
  
  # # do something...
  # N <- length(metab[1, ])
  # for (i in 1:N) {
  #   dum <- metab[, i, with = F] > -1
  #   dum <- c(dum)
  #   set_to <- unlist(metab[dum, i, with = F]/10000)
  #   set_to <- as.double(set_to)
  #   set(x = metab, i = which(dum), j = i, value = set_to)
  # }
  
  ##################################
  # Relatedness Adjustment         #
  ##################################
  
  # Final overlap of metabolites with corresponding relationship information
  overlap1 <- venn2(zz$pid,
                    ind$SampleID,
                    mytitle = "Samples for relationship adsjustment",
                    mylabels = c("Metabolite Samples \n available",
                                 "Population infos \n available"))
  
  # set the final filter for the adjustment
  filter <- (!is.element(zz$pid, ethnic_outlier))
  metab_adj <- matrix(-10000, nrow = nrow(metab[, 1, with = F]), ncol = length(metab[1, ]))
  dim(metab_adj)
  M <- length(metab[1,])
  
  # edit Carl 20181223 - restrict ind1 and ind2 to ids in adjust_me
  ind1 <- ind1[ind1$SampleID %in% pop_id, ]
  ind2 <- ind2[ind2$SampleID %in% pop_id, ]
  
  #Filter fuer Genotypen
  filt <- matrix(-1, nrow = length(ind1[,1]) + length(ind2[,1]),ncol = M)
  
  for (i in 1:M) {
    
    # i <- 1
    filter_ges <- c(filter & (metab[, i, with = F] > -10000)) # edit 21.12.18 to allow values down to -10000, e.g. normalized data
    mm <- metab[filter_ges, i, with = F]
    dat <- zz[filter_ges, ]
    vw1 <- vw[filter_ges, filter_ges]
    
    #baue Genotyp-Filter
    filt1 <- (is.element(ind1$SampleID, dat$pid))
    filt2 <- (is.element(ind2$SampleID,dat$pid) &
                !is.element(ind2$SampleID,ind1$SampleID))
    filt[,i] <- c(filt2,filt1)
    res <- lm(unlist(mm) ~ 1)$residuals
    d <- data.frame(res)
    names(d) <- "test"
    z <- polygenic(d$test, kin = vw1, d, starth2 = 0.5, quiet = TRUE)
    dum <- z$grresidual * sqrt(2)
    metab_adj[filter_ges,i] <- dum + mean(unlist(mm))
  } #i
  
  # now I should have adjusted metabolites
  metab_adjust_final <- as.data.table(metab_adj)
  setnames(metab_adjust_final, old = names(metab))
  metab_adjust_final[, pid := zz$pid]
  invisible(moveColFront(metab_adjust_final, "pid"))

  # turn -10000s back into NAs
  metab_adjust_final[, names(metab_adjust_final) := lapply(.SD, function(x){
    x[x == -10000] <- NA
    x
  })]
  
  res <- merge(x = metab_adjust_final,
               y = zz[ , .SD, .SDcols = c("pid", "FileName", "date", "SampleIndex", "SampleName", "Comment")],
               by = "pid")
  
  return(res)
} # end function