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
  t.stat  <- mediation.beta  / SE
  # t.stat2 <- mediation.beta2  / SE
  
  # when tauAdj is larger than tau we have an increased effect, which is weird and we don't want
  # therefore, invert the
  # invert.test <- ifelse(abs(tau) > abs(tauAdj), FALSE, TRUE)
  
  # when t.stat is negative, check lower tail
  # lower.tail <- ifelse(sign(t.stat) == -1, TRUE, FALSE)
  
  # calculate the p-value
  # use one-sided test and only consider positive results in the positive are of the sampling distribution
  # that means: No absolute and no multiplication by 2
  p.value <- 2*pnorm(q = abs(t.stat),
                   mean = 0,
                   sd = 1,
                   lower.tail = F)
  
  # enter into result table
  res <- data.table(
    beta.tau.tauAdj = mediation.beta,
    beta.alpha.beta = mediation.beta2,
    se = SE,
    t.statistic = t.stat,
    # t.statistic.alt = t.stat2,
    p.value = p.value
  )
  
  return(res)
}