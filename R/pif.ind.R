
#' @title PIF estimation with individual-level exposure data
#' @description Estimates the potential impact fraction (PIF) and population
#' attributable fraction (PAF) with individual-level exposure data. Only linear
#' counterfactual exposures of the form g(x) = a + b * x for the PIF are
#' supported. By default, the PAF is estimated if no counterfactual exposure
#' values are specified.
#'
#'
#' @param x vector of exposure values
#' @param beta beta coefficient in exponential relative risk
#' @param varbeta variance of beta coefficient in relative risk
#' @param a the intercept (single value) for counterfactual exposure.
#' Must be negative.
#' @param b the slope (single value) for the counterfactual exposure.
#' Must be between 0 and 1.
#' @param alpha 100*(1-alpha)\% confidence interval of PIF estimate to return
#'
#' @return A list of the following:
#' \item{\code{pif}}{point estimate of the PIF }
#' \item{\code{ci}}{100*(1-alpha)\% confidence interval of PIF estimate}
#'
#' @export
#' @import stats
#'
#' @examples
#' # Generate some exposure data
#' set.seed(1)
#' x <- rweibull(100, 1.2, 1.66)
#'
#' # Estimate PAF
#' pif.ind(x, log(1.27), 0.002)
#'
#' # Estimate PIF for a counterfactual exposure of a 1 unit decrease
#' pif.ind(x, log(1.27), 0.002, a = -1)
#'
#' # Estimate PIF for a counterfactual exposure of 25% decrease
#' pif.ind(x, log(1.27), 0.002, b = 0.75)
#'
#' @references
#' Colleen E. Chan, Rodrigo Zepeda-Tello, Dalia Camacho-Garcia-Formenti,
#' Frederick Cudhea, Rafael Meza, Eliane Rodrigues, Donna Spiegelman, Tonatiuh
#' Barrientos-Gutierrez, and Xin Zhou (2022).
#' Nonparametric Estimation of the Potential Impact Fraction and Population
#' Aggregate Fraction with Individual-Level and Aggregated Data.
#' \url{https://arxiv.org/pdf/2207.03597.pdf}
#'
pif.ind <- function(x,
                    beta,
                    varbeta,
                    a = 0,
                    b = 0,
                    alpha = 0.05){
  n <- length(x)
  if (sum(is.na(x)) != 0){
    stop("NA values supplied in x")
  }
  estpaf <- 0
  if (a != 0 & b != 0){ # both a and b are specified
    if (b > 1 | b <= 0){
      stop("b must be greater than 0 and less 1.")
    }
    if (a >= 0){
      stop("a must be a negative value or specify b value.")
    }
    gx <- a + b * x
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = ", a, " + ", b, "x",
      " and ", 100 * (1 - alpha), "% confidence interval")
    )
  } else if (a != 0){ # a is specified
    if (a >= 0){
      stop("a must be a negative value or specify b value.")
    }
    gx <- a + x
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = x - ", abs(a),
      " and ", 100 * (1 - alpha), "% confidence interval")
    )
  } else if (b != 0){ # only b is specified
    if (b >= 1 | b <= 0){
      stop("b must be between 0 and 1 or specify a nonzero intercept.")
    }
    gx <- b * x
    message(paste0(
      "Estimating PIF with counterfactual exposure g(x) = ", b, "x",
      " and ", 100 * (1 - alpha), "% confidence interval"))
  } else{
    estpaf <- 1
  }

  if (!estpaf){
    gx[gx < 0] <- 0
    if (sum(gx) == 0){
      warning("The counterfactual exposure contains all 0 values.
            The PAF will be estimated.")
      estpaf <- 1
    }
  }

  if (estpaf){
    message(paste0("Estimating PAF and ", 100 * (1 - alpha),
    "% confidence interval"))
    pif <- 1 - 1/mean(exp(beta * x))
    v.mu <- (mean( exp(beta * x)^2 ) - (mean(exp(beta * x)))^2)/n +
      varbeta * (mean(x * exp(beta * x)))^2
    v <- v.mu / (mean(exp(beta * x)))^4
  } else{
    pif <- 1 - 1/mean(exp(beta * gx))
    v.mu <- (mean( exp(beta * gx)^2 ) - (mean(exp(beta * gx)))^2)/n +
      varbeta * (mean(x * exp(beta * x)))^2
    v <- v.mu / (mean(exp(beta * x)))^4
  }
  ci <- c(pif - qnorm(1 - alpha / 2) * sqrt(v),
          pif + qnorm(1 - alpha / 2) * sqrt(v))


  return(list(pif = pif, ci = ci))
}
