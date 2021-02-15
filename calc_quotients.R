#================================================================#
# QUOTIENT CALC FUNCTION ----------------------------------------#
calc_quotients <-  function(metab, calculator = calc) {
  
  # First I create a dummy column so I don't have to worry about the "0" entries
  # metab_a1[, `:=`("0" = rep(0, .N))]
  lapply(calculator[["Quotient"]], function(Q) {
    
    # Q <- calculator[["Quotient"]][29] # Debug
    if (calculator[Quotient == (Q), Numerator2 == 0] &
        calculator[Quotient == (Q), Denominator2 == 0]) {
      toset <-
        metab[, get(calculator$Numerator1[which(calculator$Quotient == (Q))])] /
        metab[, get(calculator$Denominator1[which(calculator$Quotient == (Q))])]
    } else if (calculator[Quotient == (Q), Numerator2 != 0] &
               calculator[Quotient == (Q), Denominator2 == 0]) {
      toset <-
        (metab[, get(calculator$Numerator1[which(calculator$Quotient == (Q))])] +
           metab[, get(calculator$Numerator2[which(calculator$Quotient == (Q))])]) /
        metab[, get(calculator$Denominator1[which(calculator$Quotient == (Q))])]
    } else if (calculator[Quotient == (Q), Numerator2 == 0] &
               calculator[Quotient == (Q), Denominator2 != 0]) {
      toset <-
        metab[, get(calculator$Numerator1[which(calculator$Quotient == (Q))])] /
        (metab[, get(calculator$Denominator1[which(calculator$Quotient == (Q))])] +
           metab[, get(calculator$Denominator2[which(calculator$Quotient == (Q))])])
    } else if (calculator[Quotient == (Q), Numerator2 != 0] &
               calculator[Quotient == (Q), Denominator2 != 0]) {
      toset <-
        (metab[, get(calculator$Numerator1[which(calculator$Quotient == (Q))])] +
           metab[, get(calculator$Numerator2[which(calculator$Quotient == (Q))])]) /
        (metab[, get(calculator$Denominator1[which(calculator$Quotient == (Q))])] +
           metab[, get(calculator$Denominator2[which(calculator$Quotient == (Q))])])
    } else {
      toset <- rep("ERROR", nrow(m))
    }
    metab[, (Q) := toset]
  })
  
  # Are all all quotients present after calculation?
  if ((length(which(colnames(metab) %in% calculator$Quotient)) == length(calculator$Quotient)) == T)
    message("All Quotients have been calculated!") else message("Something seems to be missing. Please check.")
  
  # get the quotient names
  quotients <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q8", "Q9",
                 "Q11", "Q12", "Q13", "Q14", "Q15", "Q16", "Q17",
                 "Q18", "Q19", "Q20", "Q21", "Q22", "Q23", "Q24", 
                 "Q25", "Q26", "Q27", "Q28", "Q30", "Q32", "Q34",
                 "Q35",  "Q36", "Q38", "Q37")
  
  # Set the Inf quotients to the finite max
  dum1 <- metab[, lapply(.SD, max, na.rm = T), .SDcols = quotients]
  infquotients <- quotients[dum1 == Inf]
  metab[, (infquotients) := lapply(.SD, function(x) {
    
    # if (any(is.infinite(x))) {
      f <- !(is.infinite(x))
      m <- max(x[f], na.rm = T)
      x[!f] <- m
      x
    # }
  }), .SDcols = infquotients]
  
  # Turn NaN to NA
  metab[, (quotients) := lapply(.SD, function(x) {
    x[is.nan(x)] <- NA
    x}), .SDcols = quotients]
  
  # return new metab with quotients
  return(metab)
}
# END: QUOTIENT CALC FUNCTION -----------------------------------#
#================================================================#