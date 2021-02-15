cohort_specific_partial_r_squared_elimination <- function(correctionMethod,
                                                          rSquaredCutoff,
                                                          referenceData,
                                                          rawData,
                                                          allResponses,
                                                          InclHighMissings = T){
  
  # start with all covars - matching with reshaped.data in the function call chooses the right covariates for each cohort
  all.terms <- referenceData[correction == (correctionMethod) & r.squared.cutoff == (rSquaredCutoff), .(sigs = sum(sigs)), by = term]
  all.terms <- all.terms[sigs > 0, term]
  
  # record the starting terms for correct melting of all covariate columns after the loop
  melt.terms <- all.terms
  
  # seperate terms for all cohorts 
  all.terms.list <- list(a1 = all.terms,
                         b3 = all.terms,
                         sorb = all.terms)
  
  if(InclHighMissings == T){
    
    # cross-reference this list against the annotation table to exclude covars a cohort has no data for
    all.terms.list[["a1"]] <- all.terms.list[["a1"]][all.terms.list[["a1"]] %in% annot$covar.annot[a1.reason.for.exclusion %nin% "not_measured", covariate]]
    all.terms.list[["b3"]] <- all.terms.list[["b3"]][all.terms.list[["b3"]] %in% annot$covar.annot[b3.reason.for.exclusion %nin% "not_measured", covariate]]
    all.terms.list[["sorb"]] <- all.terms.list[["sorb"]][all.terms.list[["sorb"]] %in% annot$covar.annot[sorb.reason.for.exclusion %nin% "not_measured", covariate]]
    
  } else if(InclHighMissings == F) {
    
    # hm = high missings
    hm <- annot$covar.annot[a1.missing.rel >= .1 | 
                              b3.missing.rel >= .1 |
                              sorb.missing.rel >= .1,
                            covariate]
    
    # cross-reference this list against the annotation table to exclude covars a cohort has no data for
    all.terms.list[["a1"]] <- all.terms.list[["a1"]][all.terms.list[["a1"]] %nin% hm[hm != "fasting.hours"]]
    all.terms.list[["b3"]] <- all.terms.list[["b3"]][all.terms.list[["b3"]] %nin% hm[hm != "fasting.hours"]]
    all.terms.list[["sorb"]] <- all.terms.list[["sorb"]][all.terms.list[["sorb"]] %nin% hm]
    
  }
  # prime the repeat loops with all confounders
  new.confounders.list <- all.terms.list
  
  # collect cohort-based results
  my.collection <- NULL
  
  # start index for iterations
  my.iter <- 1
  
  # repeat the following code, i.e. the partial r-squared computation and elimination of covariates that do not have at least one 5% r-squared hit
  repeat{
    
    # filter for old excluded confounders
    all.terms.list <- mapply(all.terms.list,
                             new.confounders.list,
                             FUN = function(x, y){x[x %in% y]},
                             SIMPLIFY = F)
    
    # Cohort-specific partial-r-squared calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    # A1
    # filter covars based on covariates available for the cohort
    a1.terms <- all.terms.list[["a1"]]
    # a1.terms <- a1.terms[a1.terms %in% annot$covar.annot[a1.in.confounder.analysis == T, covariate]]
    
    # a1 predictors partial r squared
    a1.res <- partial_r_squared(cohort = "a1",
                                predictors = a1.terms,
                                responses = allResponses,
                                data = rawData,
                                rSquaredCutoff = rSquaredCutoff,
                                verbose = F)
    
    # B3
    # filter for old excluded confounders
    b3.terms <- all.terms.list[["b3"]]
    # b3.terms <- b3.terms[b3.terms %in% annot$covar.annot[b3.in.confounder.analysis == T, covariate]]
    
    # b3 predictors partial r squared
    b3.res <- partial_r_squared(cohort = "b3",
                                predictors = b3.terms,
                                responses = allResponses,
                                data = rawData,
                                rSquaredCutoff = rSquaredCutoff,
                                verbose = F)
    
    # SORB
    # filter for old excluded confounders
    sorb.terms <- all.terms.list[["sorb"]]
    # sorb.terms <- sorb.terms[sorb.terms %in% annot$covar.annot[sorb.in.confounder.analysis == T, covariate]]
    
    # sorb predictors partial r squared
    sorb.res <- partial_r_squared(cohort = "sorb",
                                  predictors = sorb.terms,
                                  responses = allResponses,
                                  data = rawData,
                                  rSquaredCutoff = rSquaredCutoff,
                                  verbose = F)
    
    # exclude unimportant covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    
    # define old/current confounders
    old.confounders.list <- new.confounders.list
    
    # bind together
    all.res <- rbindlist(list(a1.res, b3.res, sorb.res), use.names = T, fill = T)
    
    # special case when a covariate is not present in a cohort
    my.measure.vars <- unique(unlist(old.confounders.list))[unique(unlist(old.confounders.list)) %in% names(all.res)]
    
    # melt the covariates 
    melt.res <- melt.data.table(all.res,
                                id.vars = "cohort",
                                measure.vars = my.measure.vars,
                                value.name = "r.squared",
                                variable.name = "covariate", variable.factor = F)
    
    # plot r-squared distribution NEEDS REWORK
    # p1 <- plot_r_squared(myData = all.res,
    #                      confounders = old.confounders,
    #                      correction = correctionMethod,
    #                      rSquaredCutoff = rSquaredCutoff)
    
    # try to return plot -> wie kann ich den plot ausgeben???
    # print(p1)
    
    melt.max.res <- melt.res[ , .(max.r.squared = max(r.squared, na.rm = T)), by = .(cohort, covariate)]
    
    # change -Inf to NA
    change_in_dt(dat = melt.max.res, from = -Inf, to = NA, change_in_dat = T)
    
    # get the smalles of the maxima
    melt.max.min.res <- melt.max.res[, .SD[which.min(max.r.squared)], by = cohort]
    
    # for plot ------------------------------------------------------ 
    # record maximum of each covariate by cohort
    cohort.max.res <- all.res[, lapply(.SD, max, na.rm = T),
                              .SDcols = my.measure.vars,
                              by = .(cohort, r.squared.cutoff)]
    
    # change -Inf to NA
    change_in_dt(dat = cohort.max.res, from = -Inf, to = NA, change_in_dat = T)
    
    # add cutoff identifier
    cohort.max.res[, `:=`(correction = correctionMethod,
                          r.squared.cutoff = rSquaredCutoff,
                          iteration = my.iter)]
    
    # add to iteration
    my.iter <- my.iter + 1
    
    # add to list
    my.collection <- append(my.collection, list(cohort.max.res))
    # End: for plot -------------------------------------------------
    
    if(any(melt.max.min.res$max.r.squared < rSquaredCutoff)) {
      
      # check each list element for cutoff criteria and exclude covariate depending on criteria
      new.confounders.list <- sapply(melt.max.min.res$cohort, function(my.cohort){
        
        # get the old cohort specific confounders
        my.old <- old.confounders.list[[my.cohort]]
        
        if(melt.max.min.res[cohort == (my.cohort), max.r.squared < (rSquaredCutoff)]){
          
          # exclude covariate when r-squared is below cutoff
          my.new <- my.old[my.old %nin% melt.max.min.res[cohort == (my.cohort), covariate]]
          
        } else {
          
          # no change when r-squared not below cutoff
          my.new <- my.old  
        }
        
        # output
        return(my.new)
        
      }, USE.NAMES = T, simplify = F) # END OF APPLY
      
    } else if(all(melt.max.min.res$max.r.squared >= rSquaredCutoff)) {break} # stop!
    
    # pretty things up a little bit
    message("Finished iteration: ", my.iter, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
  } # end of repeat loop
  
  # # plot the process of the loop
  # plot.data <- rbindlist(my.collection, use.names = T, fill = T)
  # plot.data <- melt.data.table(data = plot.data,
  #                              variable.name = "covariate",
  #                              value.name = "r.squared",
  #                              id.vars = c("cohort","correction", "r.squared.cutoff", "iteration"),
  #                              measure.vars = melt.terms)
  # 
  # # draw plot
  # my.plot <- ggplot(plot.data, aes(x = factor(iteration), y = r.squared, col = covariate, shape = covariate)) + 
  #   geom_point(cex = 5) + 
  #   geom_line(aes(group = covariate)) +
  #   # scale_y_continuous(breaks = pretty_breaks(6)) +
  #   scale_y_sqrt(breaks = pretty_breaks(6)) +
  #   geom_hline(yintercept = rSquaredCutoff, col = "red", lty = "dotted") + 
  #   facet_grid(cohort~.) + 
  #   scale_shape_manual(values = unique(factor(plot.data$covariate))) + 
  #   ggtitle(label = paste0(correctionMethod, " metab~covar partial r-squared"),
  #           subtitle = paste0("Stepwise elimination of covariates with r-squared < ", rSquaredCutoff)) +
  #   xlab(label = "Iteration") +
  #   theme_carl()
  # 
  # # print for all world to see
  # print(my.plot)
  
  # max.res <- dcast.data.table(melt.max.res, formula = cohort ~ covariate, value.var = "max.r.squared")
  
  # return the old.confounders variable (no reassignment as new.confounder when max.min >= 0.05)
  return(melt.max.res)
}