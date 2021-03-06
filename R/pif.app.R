
#' @title PIF approximate estimation with mean and variance
#' @description Estimates the potential impact fraction (PIF) and population
#' attributable fraction (PAF) using only the mean and variance of the exposure
#' data. Only linear counterfactual exposures of the form g(x) = a + b * x for
#' the PIF are supported. By default, the PAF is estimated if no counterfactual
#' exposures values are specified.
#'
#' @param meanx mean of the exposure data
#' @param varx variance of the exposure data
#' @param n sample size of exposure
#' @param beta beta coefficient in exponential relative risk
#' @param varbeta variance of beta coefficient in relative risk
#' @param a the intercept (single value) for counterfactual exposure.
#' Must be negative.
#' @param b the slope (single value) for the counterfactual exposure.
#' Must be between 0 and 1.
#' @param alpha 100*(1-alpha)\% confidence interval of PAF estimate
#'
#' @return A list of the following:
#' \item{\code{pif}}{point estimate of the PIF }
#' \item{\code{ci}}{100*(1-alpha)\% confidence interval of PIF estimate}
#'
#' @export
#' @import stats
#'
#' @examples
#'
#' # Estimate PAF for an exposure distribution with mean 1.41 and variance 1.18
#' pif.app(meanx = 1.41, varx = 1.18, n = 100, beta = log(1.27),
#' varbeta = 0.002)
#'
#' # Estimate PIF for a counterfactual exposure of a 1 unit decrease
#' pif.app(meanx = 1.41, varx = 1.18, n = 100, beta = log(1.27),
#' varbeta = 0.002, a = -1)
#'
#' # Estimate PIF for a counterfactual exposure of 50% decrease
#' pif.app(meanx = 1.41, varx = 1.18, n = 100, beta = log(1.27),
#' varbeta = 0.002, b = 0.5)
#'
#' @references
#' Colleen E. Chan, Rodrigo Zepeda-Tello, Dalia Camacho-Garcia-Formenti,
#' Frederick Cudhea, Rafael Meza, Eliane Rodrigues, Donna Spiegelman, Tonatiuh
#' Barrientos-Gutierrez, and Xin Zhou (2022).
#' Nonparametric Estimation of the Potential Impact Fraction and Population
#' Aggregate Fraction with Individual-Level and Aggregated Data.
#' \url{https://arxiv.org/pdf/2207.03597.pdf}
#'
pif.app <- function(meanx,
                    varx,
                    n,
                    beta,
                    varbeta,
                    a = 0,
                    b = 0,
                    alpha = 0.05){

  estpaf <- 0
  if (a != 0 & b != 0){ # both a and b are specified
    if (b > 1 | b <= 0){
      stop("b must be greater than 0 and less 1.")
    }
    if (a >= 0){
      stop("a must be a negative value or specify b value.")
    }
    gmeanx <- a + b * meanx
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = ", a, " + ", b, "x",
      " and ", 100 * (1 - alpha), "% confidence interval")
    )
  } else if (a != 0){ # a is specified
    if (a >= 0){
      stop("a must be a negative value or specify b value.")
    }
    gmeanx <- a + meanx
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = x - ", abs(a),
      " and ", 100 * (1 - alpha), "% confidence interval")
    )
  } else if (b != 0){ # only b is specified
    if (b >= 1 | b <= 0){
      stop("b must be between 0 and 1 or specify a nonzero intercept.")
    }
    gmeanx <- b * meanx
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = ", b, "x",
      " and ", 100 * (1 - alpha), "% confidence interval"))
  } else{
    estpaf <- 1
  }

  if (!estpaf){
    if (gmeanx <= 0){
      warning("The counterfactual exposure mean is non-positive.
            The PAF will be estimated.")
      gmeanx <- 0
      estpaf <- 1
    }
  }

  if (estpaf){
    message(paste0("Estimating PAF and ", 100 * (1 - alpha),
                   "% confidence interval"))
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

