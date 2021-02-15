quotient_calc_table <- function() {
  
  #================================================================#
  # QUOTIENT CALC TABLE -------------------------------------------#
  
  # Get the data table to do the quotient calculation on
  # based on "/mnt/ifs1_projekte/genstat/02_projekte/1212_mqtl/quotienten_definitionen_2013-02-19.txt"
  q1 <- fread("/net/ifs1/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/170824_MetabAnnotTab/data/170328_quotientDefinition.csv")
  
  # I need a way to calculate the quotients based on the string in q1$Abbreviation_new
  # First get the formulas for the quotients so I know how to calculate them
  qf <- q1[["Abbreviation_new"]][grep(pattern = "Q", q1[["Abbreviation_new"]])]
  
  # Get the quotient name by splitting at the first ":" in the formula "qf"
  fp1 <- sapply(qf, function (i) substr(i, start = 1,
                                        stop = gregexpr(pattern = ":", text = qf[which(i == qf)])[[1]][1] - 1),
                USE.NAMES = F)
  
  # get the latter part of the formula
  fp2 <-
    sapply(qf, function (i)
      substr(i,start = gregexpr(pattern = ":",
                                text = qf[which(i == qf)])[[1]][1] + 1,
             stop = nchar(i)), USE.NAMES = F)
  
  # Create the basic calc data.table
  calc <- data.table(Quotient = fp1, Formula = fp2)
  
  # Try to turn the formula into a machine-readable format
  calc[, Formula := gsub(pattern = "\\/",
                         replacement = " / ",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "\\+",
                         replacement = " + ",
                         x = Formula)]
  
  # One exception due to questionable naming of metabolites
  # should have used the "idr" instead
  calc[, Formula := gsub(pattern = "OH \\+ HMG",
                         replacement = "OHHMG",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "\\(|\\)",
                         replacement = "",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "\\|",
                         replacement = "",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "C18\\:1",
                         replacement = "C181",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "C5\\:1",
                         replacement = "C51",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "AC-total",
                         replacement = "acges",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "C14\\:1",
                         replacement = "C141",
                         x = Formula)]
  calc[, Formula := gsub(pattern = "OH\\-Prol",
                         replacement = "OHProl",
                         x = Formula)]
  
  # Split up the formula so I can calculate the quotients
  calc[, `:=`(
    Numerator = tstrsplit(x = calc$Formula, split = " \\/ ")[[1]],
    Denominator = tstrsplit(x = calc$Formula, split = " \\/ ")[[2]]
  )]
  calc[, `:=`(
    Numerator1 = tstrsplit(x = calc$Numerator, split = " \\+ ")[[1]],
    Numerator2 = tstrsplit(x = calc$Numerator, split = " \\+ ")[[2]],
    Denominator1 = tstrsplit(x = calc$Denominator, split = " \\+ ")[[1]],
    Denominator2 = tstrsplit(x = calc$Denominator, split = " \\+ ")[[2]]
  )]
  
  # Replace NAs with NULLs for arithmetic processing
  calc[, names(calc) := lapply(.SD, function(x) {
    x[is.na(x)] <- 0
    x
  })]
  
  return(calc)
  # END: QUOTIENT CALC TABLE --------------------------------------#
  #================================================================#
}