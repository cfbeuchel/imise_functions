clean_sample_names <- function(dat, sample.ident, file.ident, sample.cols) {
  
  # set identifier 
  dat[, line.no := 1:.N]
  
  # get the duplicated identifiers (including NAs)
  dupli.ind <- dat[, duplis(unlist(.SD)), .SDcols = sample.ident]
  duplicates <- dat[dupli.ind, ]
  
  # set column order
  ident.positions <- grep(pattern = paste(c(sample.ident, file.ident), collapse = "|"), x = names(duplicates))
  setcolorder(x = duplicates, neworder = c(names(duplicates)[ident.positions], names(duplicates)[-ident.positions]))
  
  # Get the amount of NAs
  duplicates[, sum.NAs := sum(is.na(.SD)), .SDcols = sample.ident, by = sample.ident]
  
  # exclude NAs
  dupli.ids <- duplicates[, unlist(.SD), .SDcols = sample.ident]
  dupli.nas <- !is.na(dupli.ids)
  duplicates <- duplicates[dupli.nas, ]
  
  # get the creation date from the file name
  creation.date <- duplicates[, unique(unlist(.SD)), .SDcols = file.ident]
  creation.date.ordered <- creation.date[order(substr(creation.date, 1,8), decreasing = F)]
  
  # create ordered factor for file date
  duplicates[, file.name.ordered := factor(get(file.ident), levels = creation.date.ordered, ordered = T)]
  
  # Order, first by NAs and then by measurement date
  setorder(x = duplicates, sum.NAs, -file.name.ordered)
  
  # duplicated() only goes T for the latter element (based on index?)
  duplicates[, bad.duplicate := duplicated(get(sample.ident))]
  
  # reorder
  duplicates <- duplicates[duplis(get(sample.ident))]
  
  # rearrange output
  # set column order
  ident.positions.2 <- grep(pattern = paste(c(sample.ident, file.ident, "file.name.ordered", "sum.NAs", "bad.duplicate"), collapse = "|"), x = names(duplicates))
  setcolorder(x = duplicates, neworder = c(names(duplicates)[ident.positions.2], names(duplicates)[-ident.positions.2]))
  
  # find the lines that will be dropped given the criterium
  bad.lines <- duplicates[bad.duplicate == T, line.no]
  
  # filter the duplicated values after amount of NAs, date of survey or index
  remove.this <- dat[is.na(get(sample.ident)) | line.no %in% bad.lines, line.no]
  
  # filter all NAs and doubles
  dat <- dat[!(line.no %in% remove.this)]
  
  # return cleaned data
  return(dat)
}