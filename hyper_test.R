hyper_test <- function (
  pathway_gene,
  gezogene,
  alle_gene,
  more = T,
  unique = T)
  
{
  if (unique == T) {
    pathway_gene <- unique(pathway_gene)
    
    gezogene <- unique(gezogene)
    
    alle_gene <- unique(alle_gene)
    
  }
  
  pathway_gene_inbg = pathway_gene[pathway_gene %in% alle_gene]
  
  message(
    "Schraenke 'pathway_gene' auf 'alle_gene' ein: Using ",
    
    length(pathway_gene_inbg),
    " instead of ",
    length(pathway_gene)
  )
  
  pathway_gene = pathway_gene_inbg
  
  gezogene_inbg = gezogene[gezogene %in% alle_gene]
  
  message(
    "Schraenke 'gezogene' auf 'alle_gene' ein: Using ",
    
    length(gezogene_inbg),
    " instead of ",
    length(gezogene)
  )
  
  gezogene = gezogene_inbg
  
  if (all(pathway_gene %in% alle_gene) == F)
    
    stop("nicht alle 'pathway_gene' (alle weissen Kugeln) in 'alle_gene' (der alle_gene)")
  
  if (all(gezogene %in% alle_gene) == F)
    
    stop("nicht alle 'gezogene' (alle gezogenen Kugeln) in 'alle_gene' (der alle_gene)")
  
  if (any(is.na(pathway_gene),
          is.na(pathway_gene),
          is.na(pathway_gene)) ==
      
      T)
    
    stop("NA in den daten versteckt!")
  
  n_inPathway_inGezogen <- sum(pathway_gene %in% gezogene)
  n_inPathway_notGezogen <- length(pathway_gene) - n_inPathway_inGezogen
  n_notPathway_inGezogen <- length(gezogene) - n_inPathway_inGezogen
  n_notPathway_notGezogen <- length(alle_gene) - n_notPathway_inGezogen - n_inPathway_notGezogen - n_inPathway_inGezogen
  
  pval = stats::phyper(
    n_inPathway_inGezogen - 1,
    n_inPathway_inGezogen + n_inPathway_notGezogen,
    n_notPathway_inGezogen + n_notPathway_notGezogen,
    n_inPathway_inGezogen + n_notPathway_inGezogen,
    lower.tail = !more
  )
  
  in_gezogen <- round(
    (
      n_inPathway_inGezogen / (n_inPathway_inGezogen + n_notPathway_inGezogen)
    ) * 100, 3)
  
  in_bk <- round(
    (
      (
        n_inPathway_inGezogen + n_inPathway_notGezogen) / (
          n_inPathway_inGezogen + n_inPathway_notGezogen + n_notPathway_inGezogen + n_notPathway_notGezogen
        )
    ) * 100, 3)
  
  enr <- round(in_gezogen / in_bk, 3)
  
  mymatrix = matrix(
    c(
      n_inPathway_inGezogen,
      n_inPathway_notGezogen,
      n_notPathway_inGezogen,
      n_notPathway_notGezogen
    ),
    nrow = 2
  )
  
  or = stats::fisher.test(mymatrix)
  
  pvalfisher = or$p.value
  
  message1 = paste(
    in_gezogen,
    "% vs. ",
    in_bk,
    "% Enrichment:",
    enr,
    "OR (95%CI) =",
    signif(or$estimate, 3),
    paste0(
      "(",
      signif(or$conf.int[1], 3), "-", signif(or$conf.int[2]),
      ")"),
    sep = " "
  )
  
  message2 = paste("p hypergeomtrisch=",
                   signif(pval, 3),
                   "p fisher",
                   signif(pvalfisher, 3))
  
  message3 = paste(
    n_inPathway_inGezogen,
    "in",
    n_inPathway_inGezogen + n_notPathway_inGezogen,
    "gezogenen vs.",
    n_inPathway_inGezogen +
      n_inPathway_notGezogen,
    "in",
    n_inPathway_inGezogen + n_inPathway_notGezogen + n_notPathway_inGezogen + n_notPathway_notGezogen,
    "(grundgesamtheit)",
    sep = " "
  )
  
  message(message1)
  message(message2)
  message(message3)
  
  res = list(
    in_gezogen = in_gezogen,
    in_bk = in_bk,
    enrichment = enr,
    pval = pval,
    pval_fisher = pvalfisher,
    or = or$estimate,
    or_lower = or$conf.int[1],
    or_upper = or$conf.int[2],
    matrix = mymatrix,
    messages = c(message1, 
                 message2,
                 message3),
    compactresult = data.frame(
      in_gezogen = in_gezogen,
      in_bk = in_bk,
      enrichment = enr,
      pval = pval,
      pval_fisher = pvalfisher,
      or = or$estimate,
      or_lower = or$conf.int[1],
      or_upper = or$conf.int[2],
      Bes_Gez = mymatrix[1],
      Bes_nichtGez = mymatrix[3],
      nichtBes_Gez = mymatrix[2],
      nichtBes_nichtGez = mymatrix[4],
      matrix = paste(mymatrix, collapse = ", "),
      row.names = NULL
    )
  )
  
  res
  
}
