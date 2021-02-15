option_setup <- function(datatable_lines = 3,
                         datatable_nrow_allshow = 10,
                         datatable_colwidth = 40L
                         ){
  
  # do a garbage collection
  gc()
  
  # message("setting options(stringsAsFactors=FALSE)")
  options(stringsAsFactors = FALSE)
  
  # message("setting options(datatable.integer64= \"numeric\")")
  options(datatable.integer64 = "numeric")
  
  # message("setting options(\"repos\" = c(CRAN = \"http://cran.rstudio.com/\"))")
  options(repos = c(CRAN = "http://cran.rstudio.com/"))
  
  # message("setting options ( warnPartialMatchAttr = T ) and options ( warnPartialMatchDollar = T )")
  options(warnPartialMatchAttr = T)
  options(warnPartialMatchDollar = T)
  options(nwarnings = 1e+05)
  
  # time0 <<- Sys.time()
  # message("setting options for data.table( datatable.prettyprint.char = ", 
  #         datatable_colwidth, "\n)", "setting options for data.table( datatable.print.topn = ", 
  #         datatable_lines, "\n)", "setting options for data.table( datatable.print.nrows = ", 
  #         datatable_nrow_allshow, "\n)")
  
  options(datatable.prettyprint.char = datatable_colwidth)
  options(datatable.print.topn = datatable_lines)
  options(datatable.print.nrows = datatable_nrow_allshow)

  # Set threads to be used by data.table to something that wont crash the server
  # getDTthreads()
  setDTthreads(4)
}
