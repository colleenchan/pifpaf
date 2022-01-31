
#' @title PAF empirical
#' @description Estimates the population aggregate fraction with individual
#' level exposure data
#'
#' @param x vector of exposure values
#' @param beta beta coefficient in relative risk
#' @param varbeta variance of beta coefficient in relative risk
#' @param ci.level \% confidence interval of PAF estimate
#'
#' @return list of PAF estimate and corresponding confidence interval
#' @export
#'
#' @examples
#' x <- rweibull(100, 1.2, 1.66)
#' paf.emp(x, log(1.27), 0.002)
paf.emp <- function(x, beta, varbeta, ci.level = 0.95){
  n <- length(x)
  paf <- 1 - 1/mean(exp(beta * x))
  v.mu <- (mean( exp(beta * x)^2 ) - (mean(exp(beta * x)))^2)/n +
    varbeta * (mean(x * exp(beta * x)))^2
  v <- v.mu / (mean(exp(beta * x)))^4
  ci <- c(paf - qnorm((1 + ci.level) / 2) * sqrt(v),
          paf + qnorm((1 + ci.level) / 2) * sqrt(v))
  return(list(paf = paf, ci = ci))
}
