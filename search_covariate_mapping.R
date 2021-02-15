search_covariate_mapping <- function(mapping, all.files, cohort) {
  
  # create index to loop over
  mapping[, index := 1:.N]
  
  # search the raw data for the covariates listed in the mapping dt
  mapping.search <- sapply(mapping$index, function(i) {
    
    # i <- mapping$index[1]
    # set the confounder name and table it may be saved in
    
    if(cohort == "B3") {
      my.code <- mapping$b3.variable.code[i]
      my.table <- mapping$b3.table[i]
    } else if (cohort == "A1") {
      my.code <- mapping$a1.code[i]
      my.table <- mapping$a1.table[i]
    }
    my.confounder <- mapping$new.names[i]
    
    # check if any element is NA and return emtpy list element
    if(is.na(my.code) | is.na(my.table)) {
      
      # save an empty indicator to show that the element was not found
      final.list.element <- list(c(my.code, my.table))
      names(final.list.element) <- my.confounder
      attr(x = final.list.element[[1]], which = "NOT_FOUND") <- T
      
      # continue the search when no NA is found  
    } else {
      
      # find the target table based on the 6 digit code 
      my.table.code <- str_sub(my.table, -7, -2)
      
      # get the position of the target table
      target.file.position <- grep(pattern = my.table.code, names(all.files))
      
      # get the name of the target table 
      target.file <- names(all.files)[target.file.position]
      
      # scan the target file for the covariable-code
      target.covariate.position <- grep(pattern = my.code, x = names(all.files[[target.file]]))
      
      
      if(!identical(target.covariate.position, integer(0))) {
        
        # get the target covariate
        target.covariate <- names(all.files[[target.file]])[target.covariate.position]  
        
        # grab the columns with the identifying variables
        ident.cols <- names(all.files[[target.file]])[1:5]
        
        # grab the target column together with the identifiers
        covar.col <- all.files[[target.file]][ , .SD, .SDcols = c(na.omit(ident.cols), target.covariate)]
        
        # save the col in a list-element with the basic confounder name for later use 
        final.list.element <- list(covar.col)
        names(final.list.element) <- my.confounder
        
        # alternatively, when the column is not found, return empty element
      } else if(identical(target.covariate.position, integer(0))) {
        
        # save an empty indicator to show that the element was not found
        final.list.element <- list(c(my.code, my.table))
        names(final.list.element) <- my.confounder
        attr(x = final.list.element[[1]], which = "NOT_FOUND") <- T
      }
    }
    
    # return the resulting list elements into the result-list
    return(final.list.element)
  })
  
  return(mapping.search)
}