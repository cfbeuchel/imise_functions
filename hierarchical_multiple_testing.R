hierarchical_multiple_testing <- function(pvalues, group, adjustMethod = "BH", q1 = 0.05, q2 = 0.05){
  
  # Concept: I have pvalues of metabolite ~ covariate associations
  # This is comparable to the  Gene       ~ SNP testing scheme
  
  # 1. pvals for all covars adjusted for each metabolite
  # 2. minimum pval of each metabolite adjusted over all metabolites
  # 3. locally adjusted pval < locally adjusted minimum pvalue
  
  
  # assemble
  x <- data.table(p = pvalues, g = group)
  
  # local adjust
  x[ , local.adjust := p.adjust(p, method = adjustMethod), by = g]
  
  # get minimal locally adjusted pvalue 
  # defined as simes's method, defining the min pvalue for each group
  min.x <- x[ , .(min.p = min(local.adjust, na.rm = T)), by = g]
  
  # do global adjustment by adjusting min. pvalues for each term
  min.x[, global.adjust := p.adjust(min.p, method = adjustMethod)]
  
  # filter according to sig level q1
  min.x.filtered <- min.x[global.adjust <= q1,]
  
  # match the global pvalue to all terms
  matched.global.adjust <-  min.x.filtered[match(x$g, min.x.filtered$g), global.adjust]
  
  # add to data 
  x[ , global.adjust := matched.global.adjust]
  
  # procedure to adjust
  x[, hierarchical.fdr := local.adjust < global.adjust & global.adjust <= 0.05]
  
  return(x$hierarchical.fdr)
}
