reformate_for_matrixEQTL <- function (metab4_m_po, idvar = "pid") 
{
  print(file_preview(metab4_m_po))
  move_col_front(metab4_m_po, idvar)
  metab4_m_po <- t(metab4_m_po)
  colnames(metab4_m_po) <- metab4_m_po[idvar, ]
  metab4_m_po <- metab4_m_po[-1, ]
  myrownames <- rownames(metab4_m_po)
  metab4_m_po <- apply(metab4_m_po, 2, function(x) as.numeric(x))
  metab4_m_po <- data.table::data.table(id = myrownames, metab4_m_po)
  print(file_preview(metab4_m_po))
  metab4_m_po
}