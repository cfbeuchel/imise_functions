# covars <- c("age", "sex", "log.bmi", "diabetes.status.bi", "diabetes.status.tri",
#             "diabetes.anamnese", "diabetes.medication", "hba1c.percent",
#             "smoking.status", "fasting.8.hrs", "atc.code.g03", "cholesterol",
#             "ldl.cholesterol", "hdl.cholesterol", "white.blood.cells", "lymphocytes.percent",
#             "monocytes.percent", "hematocrit", "platelets", "reticulocytes",
#             "neutrophils.percent", "erythrocytes", "basophils.percent", "eosinophils.percent",
#             "blood.hemoglobin.level", "atc.code.c10", "fasting.hours")
# 
# cont.covars <- c("age", "log.bmi", "hba1c.percent", "cholesterol",
#                  "ldl.cholesterol", "hdl.cholesterol", "white.blood.cells",
#                  "lymphocytes.percent", "monocytes.percent", "hematocrit",
#                  "platelets", "reticulocytes", "neutrophils.percent",
#                  "erythrocytes", "basophils.percent", "eosinophils.percent",
#                  "blood.hemoglobin.level", "fasting.hours")
cov_char <- function(dat, covars, continuousCovars){
  
  # ic = is cont T/F
  ic <- covars %in% cont.covars
  
  # mean, range and iqr of continuous covars
  tmp6 <- lapply(covars[ic], function(i){
    
    tmp1 <- dat[, .(median = median(as.numeric(get(i)), na.rm = T)), by = cohort]
    tmp2 <- dat[, .(range = range(get(i), na.rm = T)), by = cohort]
    
    
    tmp3 <- dat[, .(iqr = quantile(get(i), probs = c(.25,.75),  na.rm = T)), by = cohort]
    tmp3 <- tmp3[, .(iqr = paste0(signif(iqr,3), collapse = "-")), by = cohort]
    
    # tmp3 <- dat[, .(iqr = IQR(get(i),  na.rm = T)), by = cohort]
    tmp.sd <- dat[, .(sd = sd(get(i), na.rm = T)), by = cohort]
    
    # identify min and max range
    tmp2$range.type <- rep(c("min", "max"), uniqueN(dat$cohort))
    
    # create tmp2.c[ast] for merging
    tmp2.c <- dcast(tmp2, cohort~range.type, value.var = "range")
    
    # co = cohort order
    co <- tmp1$cohort
    
    # bind together
    tmp4 <- do.call(cbind, list(tmp1[match(co, cohort), ],
                               tmp2.c[match(co, cohort), ],
                               tmp3[match(co, cohort), ],
                               tmp.sd[match(co, cohort), ]))
    
    # remove duplicated cohort columns
    tmp5 <- tmp4[, .SD, .SDcols = unique(names(tmp4))]
    
    # add covariate identifier
    tmp5[, term := (i)]
    return(tmp5)
  })
  
  # bind together
  tmp7 <- rbindlist(tmp6)
  
  # create joint range covar
  tmp7[, range := paste(round(min, 2), round(max, 2), sep = " ; ")]
  tmp7[, combined.range := paste0(round(median, 2), " [", range, "]")]
  tmp7[, combined.sd := paste0(round(median, 2), " [", round(sd, 2), "]")]
  tmp7[, combined.iqr := paste0(round(median, 2), " [", iqr, "]")]
  tmp7.c <- dcast(tmp7, term ~ cohort, value.var = c("combined.range", "combined.sd","combined.iqr", "iqr", "sd"))
  
  # get non binary variables
  non.binary <- dat[
    , sapply(.SD,
             function(i){
               
               # more than 2 levels?
               length(levels(factor(i))) > 2
             }), .SDcols = covars[!ic]
    ]
  
  # which are which
  cat.covars <- covars[!ic][non.binary]
  bin.covars <- covars[!ic][!non.binary]
  
  #=~=~=~=~=~=~=~=~=~=~=~=~=~=~=#
  # for binary covars
  tmp10 <- lapply(bin.covars, function(j){
    
    # get the total N in each cohort to calculate relative covar distribution
    # total.n <- dat[get(j), .N, by = c("cohort")]
    total.n <- dat[, .N, by = cohort]
    
    # get the levels of the factor and how often they are in the dat and put each in one thingy
    tmp8 <- dat[
      !is.na(get(j)), .N, by = c("cohort", j)
      ]
    set(tmp8, j=j, value = tmp8[,as.character(unlist(.SD)),.SDcols=j])

    # create relativ covar dists
    invisible(lapply(unique(tmp8$cohort), function(i){
      
      cohort.n <- total.n[cohort == (i), N]
      tmp8[cohort == (i), N.rel := N/(cohort.n)]
      
      }))
    
    # todo:
    tmp8 <- tmp8[
      order(
        get(j), decreasing = T
      ),
      .(ordered = get(j),
              N = paste0(round(N.rel, 3)*100, "%"),
          N.abs = N
      ),
      by = cohort
      ][ordered == unique(ordered)[1]
      # ][
        ,
        .(   all.n = paste0(N.abs, " (", N, ")"),
          all.lvls = paste(ordered, collapse = ", ")
        ),
        by = .(cohort)
        ]
    tmp8[, term := j]
    
    # tmp8 <- tmp8[
    #     order(
    #       get(j), decreasing = T
    #       ),
    #     .(ordered = get(j),
    #             N = paste0(round(N.rel, 3)*100, "%")
    #       ),
    #     by = cohort
    #     ][
    #       ,
    #       .(   all.n = paste(N, collapse = " / "),
    #         all.lvls = paste(ordered, collapse = ", ")
    #         ),
    #       by = cohort
    #       ]
    # tmp8[, term := j]

    return(tmp8)
  })
  tmp11 <- rbindlist(tmp10)
  tmp11[, type := "discrete"]
  tmp11.c <- dcast(tmp11, term ~ cohort, value.var = c("all.lvls", "all.n"))
  
  #=~=~=~=~=~=~=~=~=~=~=~=~=~=~=#
  # for categorical covars
  tmp12 <- lapply(cat.covars, function(j){
    
    # get the total N in each cohort to calculate relative covar distribution
    # total.n <- dat[!is.na(get(j)), .N, by = c("cohort")]
    total.n <- dat[, .N, by = cohort]
    
    # get the levels of the factor and how often they are in the dat and put each in one thingy
    tmp8 <- dat[
      !is.na(get(j)), .N, by = c("cohort", j)
      ]
    
    # create relativ covar dists
    invisible(lapply(unique(tmp8$cohort), function(i){
      
      cohort.n <- total.n[cohort == (i), N]
      tmp8[cohort == (i), N.rel := N/(cohort.n)]
      
    }))
    
    tmp8 <- tmp8[
        order(
          get(j), decreasing = T
          ),
        .(ordered = get(j),
            N.abs = N,
            N.rel = paste0(round(N.rel, 3)*100, "%")
          ),
        by = cohort
        ][
          , .(cohort,
              ordered,
              N.abs,
              N.rel,
              N.both = paste0(N.abs, " (", N.rel, ")"))
        ][
          ,
          .(   all.n = paste(N.both, collapse = " / "),
            all.lvls = paste(ordered, collapse = ", "),
                expl = c("current, ever, never")
            ),
          by = cohort
          ]
    
    tmp8[, term := j]
    
    return(tmp8)
  })
  tmp13 <- rbindlist(tmp12)
  tmp13.c <- NA
  if(nrow(tmp13)!=0){
    tmp13[, type := "categorical"]
    tmp13.c <- dcast(tmp13, term ~ cohort, value.var = c("all.lvls", "all.n"))
  }
  
  # bind together
  res <- list(continuous = tmp7.c,
              discrete = tmp11.c,
              categorical = tmp13.c)
  return(res)
}







