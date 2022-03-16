
#' @title PIF estimation with individual-level data
#' @description Estimates the potential impact fraction and population
#' aggregate fraction with individual-level exposure data
#'
#' @param x vector of exposure values
#' @param beta beta coefficient in relative risk
#' @param varbeta variance of beta coefficient in relative risk
#' @param estpaf boolean on whether to compute PIF
#' @param a the intercept (single value) for counterfactual exposure.
#' Must be negative.
#' @param b the slope (single value) for the counterfactual exposure.
#' Must be between 0 and 1.
#' @param alpha 100*(1-alpha)\% confidence interval of PIF estimate
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
#' pif.ind(x, log(1.27), 0.002, estpaf = TRUE)
#'
#' # Estimate PIF for a counterfactual exposure of a 1 unit decrease
#' pif.ind(x, log(1.27), 0.002, a = -1)
#'
#' # Estimate PIF for a counterfactual exposure of 50% decrease
#' pif.ind(x, log(1.27), 0.002, b = 0.5)
pif.ind <- function(x,
                    beta,
                    varbeta,
                    estpaf = FALSE,
                    a = 0,
                    b = 1,
                    alpha = 0.05){
  n <- length(x)

  # argument checks
  if (a == 0 & b == 1 & !estpaf){
     stop("Must specify estpaf = TRUE or supply value(s) to a and/or b to
          compute PIF.")
  }
  if ((a != 0 | b != 1) & estpaf){
     stop("a and/or b value specified but estpaf = TRUE.")
  }

  if (!estpaf){
    if (a != 0 & b != 1){ # both a and b are specified
      if (b >= 1 | b <= 0){
        stop("b must be greater than 0 and less 1.")
      }
      message(paste0(
        "Estimating PIF with counterfactual exposure g(x) = ", a, " + ", b, "x",
        " and ", 100 * (1 - alpha), "% confidence interval")
        )
    } else if (a != 0){ # a is specified
      if (a >= 0){
        stop("a must be a negative value or specificy b value.")
      }
      message(paste0(
        "Estimating PIF with counterfactual exposure g(x) = x - ", abs(a),
        " and ", 100 * (1 - alpha), "% confidence interval")
      )
    } else{ # only b is specified
      if (b >= 1 | b <= 0){
        stop("b must be between 0 and 1 or specify a nonzero intercept.")
      }
      message(paste0(
        "Estimating PIF with counterfactual exposure g(x) = ", b, "x",
        " and ", 100 * (1 - alpha), "% confidence interval"))
    }
  }

  gx <- a + b * x
  gx[gx < 0] <- 0
  if (sum(gx) == 0){
    warning("The counterfactual exposure contains all 0 values.
            The PAF will be estimated.")
    estpaf <- 1
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
