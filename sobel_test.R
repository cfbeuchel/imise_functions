sobel_test <- function(tau, # beta_X Y~X
                       tauAdj, # beta_X Y~X+M
                       alpha, # beta_X M~X
                       seAlpha, # SE_X M~X
                       beta, # beta_M Y~X+M
                       seBeta # SE_M Y~X+M 
) {
  
  # first calculate the SE
  SE <- sqrt( alpha^2 * seBeta^2 + beta^2 * seAlpha^2 )
  
  # calculate the Standard Error term
  mediation.beta <- tau - tauAdj
  mediation.beta2 <- alpha * beta
  
  # if(signif(mediation.beta,4)!=signif(mediation.beta2,4)){
  #   message("Mediation effect calcuated via tau - tau' is unequal to calculation via alpha*beta!")
  #   }
  
  # calculate the test statistic
  # t.stat  <- mediation.beta  / SE
  t.stat <- mediation.beta2  / SE
  
  # when tauAdj is larger than tau we have an increased effect, which is weird and we don't want
  # therefore, invert the test
  # edit 200714: obsolete due to changed calculation of total effect!
  
  
  # total effect is always larger than direct effect
  # whe you calculate total effect as alpha*beta + tauAdj
  # invert.test <- ifelse(abs(tau) > abs(tauAdj), FALSE, TRUE)
  
  # 200828 - invert only when sign in not in line
  # invert.test <- ifelse(sign(tau) == sign(tauAdj), FALSE, TRUE)
  # invert.test <- ifelse(sign(alpha*beta+tauAdj) == sign(tauAdj), FALSE, TRUE)
  invert.test <- (sign(alpha*beta) == sign(tauAdj)) == FALSE
  
  # should I add something for different effect signs?
  
  # when t.stat is negative, check lower tail
  lower.tail <- ifelse(sign(t.stat) == -1, TRUE, FALSE)
  
  # calculate the p-value
  # use one-sided test and only consider positive results in the positive are of the sampling distribution
  # that means: No absolute and no multiplication by 2
  p.value <- pnorm(q = t.stat,
                   mean = 0,
                   sd = 1,
                   # lower.tail = FALSE)
                   lower.tail = ifelse(invert.test, !lower.tail, lower.tail))
  
  # calculate assymetric CI
  med.ci <- RMediation::medci(mu.x = alpha,
                              mu.y = beta,
                              se.x = seAlpha,
                              se.y = seBeta,
                              alpha = 0.05,
                              plot = F,
                              plotCI = F, 
                              type = "dop")
  
  
  # get the product t-statistic
  t.stat.prod <- (alpha/seAlpha) * (beta/seBeta)
  
  prod.p <- pprodnormal(
    # q = 0,
    # q = abs((alpha/seAlpha) * (beta/seBeta)),
    q = t.stat.prod,
    mu.x = 0,
    mu.y = 0,
    se.x = 1,
    se.y = 1,
    # mu.x = alpha,
    # se.x = seAlpha,
    # mu.y = beta,
    # se.y = seBeta, 
    lower.tail = ifelse(invert.test, !lower.tail, lower.tail),
    type = "dop"
  )
  
  prod.p <- prod.p$p
  
  
  # enter into result table
  res <- data.table(
    beta.tau.tauAdj = mediation.beta,
    beta.alpha.beta = mediation.beta2,
    se = SE,
    t.statistic = t.stat,
    p.value = p.value,
    ci.95.lower = med.ci$`97.5% CI`[1],
    ci.95.upper = med.ci$`97.5% CI`[2],
    t.statistic.prod = t.stat.prod,
    p.value.prod = prod.p
  )
  
  return(res)
}
