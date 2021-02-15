all_annotation_tables <- function(mount.point = "/net/ifs1/san_projekte/projekte/") {
  # dt with locations and file names
  data.loc <- data.table::data.table(annot.table = c("metab", "covar", "gx", "ind"),
                                     annot.dir = c(paste0(mount.point, "02_projekte/1703_ge_metab_a1_b3_sorbs/170824_MetabAnnotTab/results/"),
                                                   paste0(mount.point, "02_projekte/1703_ge_metab_a1_b3_sorbs/171113_CovarAnnotTab/results/"),
                                                   paste0(mount.point, "02_projekte/1703_ge_metab_a1_b3_sorbs/171121_GxProbeAnnotTab/results/"),
                                                   paste0(mount.point, "02_projekte/1703_ge_metab_a1_b3_sorbs/171107_IndAnnotTab/results/")),
                                     file.name = c("_a1_b3_sorb_MetaboliteAnnotation.csv",
                                                   "_a1_b3_sorb_CovariateAnnotation.csv",
                                                   "_a1_b3_sorb_GxProbeAnnotation.csv",
                                                   "a1_b3_sorb_IndividualAnnotation.csv"))
  
  # find newest annotation file based on annot.dor and file.name
  for (i in data.loc$annot.table) {
    data.table::set(x = data.loc,
                    i = which(data.loc$annot.table == i),
                    j = "newest.annot.file", value = newest_file(look_for = data.loc[annot.table == i, file.name],
                                                                 directory = data.loc[annot.table == i, annot.dir]))
  }
  
  # load all annotation tables into list
  all.tables <- list(
    # metab annotation
    metab.annot = data.loc[annot.table == "metab", data.table::fread(paste0(annot.dir, newest.annot.file))],
    
    # covar annotation
    covar.annot = data.loc[annot.table == "covar", data.table::fread(paste0(annot.dir, newest.annot.file))],
    
    # gx annotation
    gx.annot = data.loc[annot.table == "gx", data.table::fread(paste0(annot.dir, newest.annot.file))],
    
    # individual annotation
    ind.annot = data.loc[annot.table == "ind", data.table::fread(paste0(annot.dir, newest.annot.file))]
  )
  
  # return list
  return(all.tables)
}
