
#' @title PIF approximate estimation
#' @description Estimates the potential impact fraction (PIF) and population
#' aggregate fraction (PAF) using only the mean and variance of the exposure
#' data. By default, the PAF is estimated if no counterfactual exposures values
#'are specified.
#'
#' @param meanx mean of the exposure data
#' @param varx variance of the exposure data
#' @param n sample size of exposure
#' @param beta beta coefficient in relative risk
#' @param varbeta variance of beta coefficient in relative risk
#' @param a the intercept (single value) for counterfactual exposure.
#' Must be negative.
#' @param b the slope (single value) for the counterfactual exposure.
#' Must be between 0 and 1.
#' @param alpha 100*(1-alpha)\% confidence interval of PAF estimate
#'
#' @return list of PIF estimate and corresponding confidence interval
#' @export
#' @import stats
#'
#' @examples
#' # Generate exposure data
#' set.seed(1)
#' x <- rweibull(100, 1.2, 1.66)
#'
#' # Estimate PAF
#' pif.app(mean(x), var(x), length(x), log(1.27), 0.002)
#'
#' # Estimate PIF for a counterfactual exposure of a 1 unit decrease
#' pif.app(mean(x), var(x), length(x), log(1.27), 0.002, a = -1)
#'
#' # Estimate PIF for a counterfactual exposure of 50% decrease
#' pif.app(mean(x), var(x), length(x), log(1.27), 0.002, b = 0.5)
pif.app <- function(meanx,
                    varx,
                    n,
                    beta,
                    varbeta,
                    a = 0,
                    b = 1,
                    alpha = 0.05){

  if (a != 0 & b != 1){ # both a and b are specified
    if (b >= 1 | b <= 0){
      stop("b must be greater than 0 and less 1.")
    }
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = ", a, " + ", b, "x",
      " and ", 100 * (1 - alpha), "% confidence interval")
    )
    estpaf <- 0
  } else if (a != 0){ # a is specified
    if (a >= 0){
      stop("a must be a negative value or specificy b value.")
    }
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = x - ", abs(a),
      " and ", 100 * (1 - alpha), "% confidence interval")
    )
    estpaf <- 0
  } else if (b != 1){ # only b is specified
    if (b >= 1 | b <= 0){
      stop("b must be between 0 and 1 or specify a nonzero intercept.")
    }
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = ", b, "x",
      " and ", 100 * (1 - alpha), "% confidence interval"))
    estpaf <- 0
  } else{
    estpaf <- 1
  }

  gmeanx <- a + b * meanx
  if (gmeanx <= 0){
    warning("The counterfactual exposure mean is non-positive.
            The PAF will be estimated.")
    gmeanx <- 0
    estpaf <- 1
  }

  if (estpaf){
    pif <- 1 -  1 / (exp(beta * meanx) * (1 + 0.5 * sqrt(varx) * beta^2))
    grad <- c(beta / (exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varx))),
              beta^2 * exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varx))^2 /
                (4 * sqrt(varx)),
              (meanx + beta * (1 + 0.5 * beta * meanx) * sqrt(varx)) /
                (exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varx))^2)
    )
  } else{
    pif <- 1 - exp(beta * gmeanx) * (1 + 0.5 * sqrt(varx) * beta^2 * b^2) /
      (exp(beta * meanx) * (1 + 0.5 * sqrt(varx) * beta^2))

    grad <- c(-(exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varx)) *
                  (beta * b * exp(gmeanx) +
                     0.5 * beta^3 * b^3 * sqrt(varx) * exp(beta * gmeanx)) -
                  exp(beta * gmeanx) * (1 + 0.5 * sqrt(varx)) * beta^2 * b^2 *
                  (beta * exp(beta * meanx)) *
                  (1 + 0.5 * beta^2 * sqrt(varx))) /
                (exp(beta * meanx)^2 * (1 + 0.5 * beta^2 * sqrt(varx))^2),
               -(exp(beta * meanx) * (1 + 0.5 * beta * sqrt(varx)) * 0.25 *
                exp(beta * gmeanx) * varx^(-1/2) * beta^2 * b^2 -
                exp(beta * gmeanx) * (1 + 0.5 * sqrt(varx) * beta^2 * b^2) *
                exp(beta * meanx) * 0.25 * beta^2 * varx^(-1/2)) /
                (exp(beta * meanx)^2 * (1 + 0.5 * beta^2 * sqrt(varx))^2),
              -(exp(beta * meanx) * (1 + 0.5 * beta^2 * sqrt(varx)) *
                  (exp(beta * gmeanx) *
                     (gmeanx + sqrt(varx) * beta * b^2 +
                        0.5 * beta^2 * b^2 * gmeanx * sqrt(varx))) -
                  (exp(beta * gmeanx) * (1 + 0.5 * sqrt(varx)) * beta^2 * b^2) *
                  (exp(beta * meanx) * beta * sqrt(varx) +
                     (1 + 0.5 * beta^2 * sqrt(varx) * meanx *
                        exp(beta * meanx)))) /
                (exp(beta * meanx)^2 * (1 + 0.5 * beta^2 * sqrt(varx))^2))
  }
  Sigma <- diag(c(varx/n,
                  (3*varx^2)/n - (n-3)*varx^(3/2)/(n*(n-1)),
                  varbeta))
  v <- grad %*% Sigma %*% grad
  ci <- c(pif - qnorm(1 - alpha / 2) * sqrt(v),
          pif + qnorm(1 - alpha / 2) * sqrt(v))
  return(list(pif = pif, ci = ci))
}

