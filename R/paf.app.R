
#' @title PAF approximate
#' @description Estimates the population aggregate fraction with the mean and
#' variance of the exposure data
#'
#' @param meanx mean of the exposure data
#' @param varx variance of the exposure data
#' @param n sample size of exposure
#' @param beta beta coefficient in relative risk
#' @param varbeta variance of beta coefficient in relative risk
#' @param ci.level \% confidence interval of PAF estimate
#'
#' @return list of PAF estimate and corresponding confidence interval
#' @export
#'
#' @examples
#' x <- rweibull(100, 1.2, 1.66)
#' paf.app(mean(x), var(x), length(x), log(1.27), 0.002)
paf.app <- function(meanx, varx, n, beta, varbeta, ci.level = 0.95){
  paf <- 1 -  1/ (exp(beta * meanx) * (1 + 0.5 * sqrt(varx) * beta^2))
  grad <- c(beta / (exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varbeta))),
            beta^2 * exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varx))^2 /
              (4 * sqrt(varx)),
            (meanx + beta * (1 + 0.5 * beta * meanx) * sqrt(varx)) /
              (exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varx))^2)
  )
  Sigma <- diag(c(varx/n,
                  (3*varx^2)/n - (n-3)*varx^(3/2)/(n*(n-1)),
                  varbeta))
  v <- grad %*% Sigma %*% grad
  ci <- c(paf - qnorm((1 + ci.level) / 2) * sqrt(v),
          paf + qnorm((1 + ci.level) / 2) * sqrt(v))
  return(list(paf = paf, ci = ci))
}

