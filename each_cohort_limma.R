each_cohort_limma <- function(
  gx.dat,
  gx.annot,
  metab.dat,
  metabolites,
  cohortID,
  covariates = c()
){
  
  # merge to eset
  for(i in metabolites) {
    pData(gx.dat)[,i] <- metab.dat[match_hk(rownames(pData(gx.dat)), metab.dat$id),get(i)]
  }
  
  # check class
  showNA(pData(gx.dat))
  showClassDF(pData(gx.dat)) 
  
  # calculate assoc
  gx.res <- calculateLimma(todovars = metabolites,
                           conf = covariates,
                           eset = gx.dat,
                           doscale = F,
                           doconfint = T,
                           save_intermediate = F
  )
  
  # calculate fdr stats
  ge_assoc_tables <- lapply(
    gx.res,
    function(listeintrag) annotiereFDReGx(
      ge_listeneintrag = listeintrag,
      gx_probeanno = gx.annot,
      ilmn_colname = "ilmn",
      entrez_colname = "EntrezReannotated",
      hgnc_colname = "SymbolReannotated_orgHsEg")
  )
  
  # summarise
  ge_assoc_summaries <- lapply(
    ge_assoc_tables,
    function(ge_table) makeFDRstatsGenelevel(
      pval = ge_table$P.Value,
      genes = ge_table$EntrezReannotated,
      logFC = ge_table$logFC,
      qval_cutoffs = c(0.05, 0.2),
      showUpDownSeparately = T,
      showVennUpDownBoth = T
    )
  )
  
  # bind together
  ge_assoc_summaries <- rbindlist(ge_assoc_summaries)
  ge_assoc_summaries$metabolite <- names(ge_assoc_tables)
  res <- rbindlist(ge_assoc_tables)
  res[,cohort := (cohortID)]
  ge_assoc_summaries[,cohort := (cohortID)]
  
  # create hierarchical fdr stats
  res.hfdr <- res[, addHierarchFDR(
    pvalues = P.Value,
    categs = var,
    fdrmethod_level1 = "BH",
    fdrmethod_level2 = "BH",
    correctionLevel1 = "BB")]
  
  # add the stats to the results
  res$hfdr.1 <- res.hfdr$fdr_level1
  res$hfdr.2 <- res.hfdr$fdr_level2
  res$sig.hfdr <- res.hfdr$hierarch_fdr5proz
  
  return(list(
    res = res,
    summary = ge_assoc_summaries
  ))
}
