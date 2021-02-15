define_categorical_metabolites <- function(metabolites, to.transform, no.of.categories = 20, plots = F) {
  
  new.categs <- invisible(sapply(to.transform, function(mymetab) {
    
    # mymetab <- to.transform[1]
    # name of category
    mycateg <- paste0(mymetab, ".categorical")
    
    # create dummy dt for calculations
    categ.dummy <- data.table::data.table(original = metabolites[[mymetab]])
    
    # do not create a seperate
    categ.dummy[ , new.categs := cut(original, breaks = seq(min(original, na.rm = T), max(original, na.rm = T), length.out = no.of.categories + 1), include.lowest = T)]
    categ.dummy[ , ordered.categs := ordered(new.categs)]
    
    #================================================================
    # (OLD, WITH SEPERATE ZERO CATEGORY) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # create categories for values > 0 based on range 
    # categ.dummy[original > 0, new.categs := cut(original, breaks = seq(0, max(original), length.out = no.of.categories + 1), include.lowest = T)]
    # categs.created <- categ.dummy[, levels(unlist(new.categs))]
    # categ.dummy[original == 0, new.categs := factor(0, labels =  "[0]")]
    # categ.dummy[ , ordered.categs := ordered(new.categs, levels = c("[0]", (categs.created)))]
    #================================================================
    
    # I want categories to have at least 50 entries and will therefore merge them until all have the desired number of entries
    no.of.small.categs <- categ.dummy[, .N, by = ordered.categs][N <= 50, length(ordered.categs)]
    
    while(no.of.small.categs > 0) {
      
      # do conditional merging when one category contains < 50 values
      # also skip any NA category. NAs remain NAs
      small.categs <- categ.dummy[!is.na(ordered.categs), .N, by = ordered.categs][N <= 50, as.character(na.omit(ordered.categs))]
      
      # also add empty categories to small categories
      empty.categs <- levels(categ.dummy$ordered.categs)[levels(categ.dummy$ordered.categs) %nin% unique(categ.dummy$ordered.categs)]
      
      # create dummy to calculate on 
      small.dummy <- data.table(small.categs = c(empty.categs, small.categs))
      
      # find the highest categ index
      highest.categ <- length(levels(categ.dummy$ordered.categs))
      
      # find the categories following the small categories, which will be merged
      categs.to.merge.to <- match(small.dummy$small.categs,  levels(categ.dummy$ordered.categs)) + 1
      
      # check if the last category wants to merge to the non-existing next category
      if(any(categs.to.merge.to > highest.categ)) {
        
        # change the category the last category merges to from the next to the previous
        last.categ.merge.to <- categs.to.merge.to[categs.to.merge.to > highest.categ]
        
        # substract 2 to get from the next to the previous category
        last.categ.merge.to.previous <- last.categ.merge.to - 2
        
        # change entry 
        categs.to.merge.to[categs.to.merge.to > highest.categ] <- last.categ.merge.to.previous
      }
      
      # enter category to merge to to small.dummy 
      small.dummy[ , merge.to := levels(categ.dummy$ordered.categs)[categs.to.merge.to]]
      
      # check for categories that want to merge that already will have been merged to another category and skip them
      small.dummy[ , skip := small.categs %in% merge.to]
      
      # build in special case when two categs want to join each other and block the merging
      if(all(small.dummy$skip == T)) {
        small.dummy[, skip := c(F,T)]
      }
      
      # for as long as the small.dummy dt is larger than nrow = 1, this will work
      if(nrow(small.dummy) > 1) {
        
        # create new category
        new.categ <- paste0(small.dummy$small.categs, small.dummy$merge.to)
        new.categ <- gsub(pattern = ",(.*),", replacement = ",", new.categ)
        small.dummy[ , new.categ := new.categ]
        
        # special case of last merge
      }  else {
        
        # determine the smaller category
        new.categ <- paste0(small.dummy$merge.to, small.dummy$small.categs)
        new.categ <- gsub(pattern = ",(.*),", replacement = ",", new.categ)
        small.dummy[ , new.categ := new.categ]
      }
      
      # try to rename the levels of the factor and thus change the entries in the dt
      my.levels <- levels(categ.dummy$ordered.categs)
      
      # match the new categories over the the original levels
      # write over both categories to be replaced. The one that was designated to be merged and the one it was supposed to be merged to
      my.levels[match(small.dummy[skip == F, small.categs], my.levels)] <- small.dummy[skip == F, new.categ]
      my.levels[match(small.dummy[skip == F, merge.to], my.levels)] <- small.dummy[skip == F, new.categ]
      
      # now reassign the new levels to the factor
      levels(categ.dummy$ordered.categs) <- my.levels
      
      # check the no of categories again
      new.no.of.small.categs <- categ.dummy[!is.na(ordered.categs), .N, by = ordered.categs][N <= 50, length(ordered.categs)]
      no.of.small.categs <- new.no.of.small.categs
    }
    
    if(plots == T) {
      
      # plotting parameters
      par(mfrow = c(1,2))
      
      # plot 10 categs
      hist(categ.dummy$original, breaks = 50, main = "New categories", xlab = mymetab)
      categ.margins <- seq(0, max(categ.dummy$original, na.rm = T), length.out = no.of.categories + 1)
      abline( v = categ.margins , lwd = 2, col = "red")
      
      # plot new categories
      hist(categ.dummy$original, breaks = 50, main = "Adjusted categories", xlab = mymetab)
      adj.categ.margins <- categ.dummy[, range(original), by = ordered.categs][ , unique(na.omit(V1))]
      abline( v = adj.categ.margins , lwd = 2, col = "red")
      
    }
    
    # return categ.dummy with correct name
    res <- list(categ.dummy$ordered.categs)
    names(res) <- mycateg
    return(res)
  }, USE.NAMES = F))
  
  # reset plotting parameters
  if(all(par()$mfrow != c(1,1))) {
    par(mfrow = c(1,1))
  }
  
  # turn into dt
  setDT(new.categs)
  
  # add categorical metabolites to metabolite data
  metabolites[ , names(new.categs) := new.categs]
}
