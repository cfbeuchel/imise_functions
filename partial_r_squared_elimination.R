partial_r_squared_elimination <- function(correctionMethod, rSquaredCutoff, referenceData, rawData){
  
  # start with all covars - matching with reshaped.data in the function call chooses the right covariates for each cohort
  all.terms <- referenceData[correction == (correctionMethod) & r.squared.cutoff == (rSquaredCutoff), .(sigs = sum(sigs)), by = term]
  all.terms <- all.terms[sigs > 0, term]
  
  # record the starting terms for correct melting of all covariate columns after the loop
  melt.terms <- all.terms
  
  # prime the repeat loops with all confounders
  new.confounders <- all.terms
  
  # collect cohort-based results
  my.collection <- NULL
  
  # start index for iterations
  my.iter <- 1
  
  # repeat the following code, i.e. the partial r-squared computation and elimination of covariates that do not have at least one 5% r-squared hit
  repeat{
    
    # filter for old excluded confounders
    all.terms <- all.terms[all.terms %in% new.confounders]
    
    # Cohort-specific partial-r-squared calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    # A1
    # filter covars based on covariates available for the cohort
    a1.terms <- all.terms[all.terms %in% annot$covar.annot[a1.in.confounder.analysis == T, covariate]]
    
    # a1 predictors partial r squared
    a1.res <- partial_r_squared(cohort = "a1",
                                predictors = a1.terms,
                                responses = unique(my.data$metab),
                                data = rawData,
                                rSquaredCutoff = rSquaredCutoff,
                                verbose = F)
    
    # B3
    # filter for old excluded confounders
    b3.terms <- all.terms[all.terms %in% annot$covar.annot[b3.in.confounder.analysis == T, covariate]]
    
    # b3 predictors partial r squared
    b3.res <- partial_r_squared(cohort = "b3",
                                predictors = b3.terms,
                                responses = unique(my.data$metab),
                                data = rawData,
                                rSquaredCutoff = rSquaredCutoff,
                                verbose = F)
    
    # SORB
    # filter for old excluded confounders
    sorb.terms <- all.terms[all.terms %in% annot$covar.annot[sorb.in.confounder.analysis == T, covariate]]
    
    # sorb predictors partial r squared
    sorb.res <- partial_r_squared(cohort = "sorb",
                                  predictors = sorb.terms,
                                  responses = unique(my.data$metab),
                                  data = rawData,
                                  rSquaredCutoff = rSquaredCutoff,
                                  verbose = F)
    
    # exclude unimportant covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    
    # define old/current confounders
    old.confounders <- new.confounders
    
    # bind together
    all.res <- rbindlist(list(a1.res, b3.res, sorb.res), use.names = T, fill = T)
    
    # plot r-squared distribution
    p1 <- plot_r_squared(myData = all.res,
                         confounders = old.confounders,
                         correction = correctionMethod,
                         rSquaredCutoff = rSquaredCutoff)
    
    # try to return plot -> wie kann ich den plot ausgeben???
    print(p1)
    
    # record maximum of each covariate by cohort
    cohort.max.res <- all.res[, lapply(.SD, max, na.rm = T), .SDcols = old.confounders, by = .(cohort, r.squared.cutoff)]
    
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
    
    # output the maximum of each covariate
    max.res <- all.res[, unlist(lapply(.SD, max, na.rm = T)), .SDcols = old.confounders]
    
    # which is the smallest of the maximum?
    max.min <- max.res[which(max.res == min(max.res))]
    
    # verbose maximum explained variance per covariate
    message("Maximum explained variance for each covariate were:")
    message(paste0(names(max.res),collape = ", "))
    message(paste(round(max.res,3), collapse = "      "))
    
    # verbose last excluded confounder with amount of variance explained
    message(paste0("Min confounder is ", names(max.min), " with ", round(max.min * 100, 2), "% explained variance."))
    
    # remove confounder from list when r2 is smaller than 5%
    if(max.min < rSquaredCutoff) {
      
      # set the new confounder set without min.max
      new.confounders <- old.confounders[old.confounders %nin% names(max.min)]
      
      # verbose information about the old/new confounders the loop is repeated with
      message(paste0(length(old.confounders), " old confounders were ", paste(old.confounders, collapse = ", ")))
      message(paste0(length(new.confounders), " new confounders are ", paste(new.confounders, collapse = ", ")))
      
      # when max.min >= 0.05
    } else {
      
      # different message in case max.min is >0.05
      message(paste0("New confounders = old confounders! Loop is done!\nTotal number is: ", length(old.confounders)))
      message(paste0("Final confounders after exclusion are ", paste0(old.confounders, collapse = ", ")))
    } # end of if/else statement
    
    # stop loop when all covariates explain 6% of response variance in at least 1 case 
    if(max.min >= rSquaredCutoff){break} # !!! break statement
    
    # pretty things up a little bit
    message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
  } # end of repeat loop
  
  # plot the process of the loop
  plot.data <- rbindlist(my.collection, use.names = T, fill = T)
  plot.data <- melt.data.table(data = plot.data,
                               variable.name = "covariate",
                               value.name = "r.squared",
                               id.vars = c("cohort","correction", "r.squared.cutoff", "iteration"),
                               measure.vars = melt.terms)
  
  # draw plot
  my.plot <- ggplot(plot.data, aes(x = factor(iteration), y = r.squared, col = covariate, shape = covariate)) + 
    geom_point(cex = 5) + 
    scale_y_continuous(breaks = pretty_breaks(6)) +
    geom_hline(yintercept = rSquaredCutoff, col = "red", lty = "dotted") + 
    facet_grid(cohort~.) + 
    scale_shape_manual(values = unique(factor(plot.data$covariate))) + 
    ggtitle(label = paste0(correctionMethod, " metab~covar partial r-squared"),
            subtitle = paste0("Stepwise elimination of covariates with r-squared < ", rSquaredCutoff)) +
    xlab(label = "elimination round") +
    theme_carl()
  
  # print for all world to see
  print(my.plot)
  
  
  # return the old.confounders variable (no reassignment as new.confounder when max.min >= 0.05)
  return(max.res)
}