
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pifpaf

<!-- badges: start -->
<!-- badges: end -->

The goal of pifpaf is to estimate potential impact fractions and
population attributable fractions using individual-level or aggregate
data.

## Installation

You can install the most updated version of pifpaf as follows:

``` r
library(devtools)
devtools::install_github("colleenchan/pifpaf")
```

## Example

This is a basic example which shows you how to estimate the PAF when
individual-level exposure data is available to the user, where the
relative risk function takes the form RR(*x*) = exp (*x*) = 1.27 with
variance 0.002.

``` r
library(pifpaf)
# generate some exposure data
set.seed(1)
x <- rweibull(1000, 1.2, 1.66) 

# Estimate PAF
pif.ind(x, beta = log(1.27), varbeta = 0.002, estpaf = TRUE) 
#> Computing PAF and 95% confidence interval
#> $pif
#> [1] 0.3447482
#> 
#> $ci
#> [1] 0.2273862 0.4621102
```

For the PIF, we allow counterfactual exposures of the form
*g*(*x*) = *a* + *b**x*.

``` r
# Estimate PIF for a counterfactual exposure of a 1 unit decrease
pif.ind(x, beta = log(1.27), varbeta = 0.002, a = -1)
#> Computing PIF and 95% confidence interval counterfactual exposure g(x) = x-1
#> $pif
#> [1] 0.198881
#> 
#> $ci
#> [1] 0.08201511 0.31574683

# Estimate PIF for a counterfactual exposure of 50% decrease
pif.ind(x, beta = log(1.27), varbeta = 0.002, b = 0.5)
#> Computing PIF and 95% confidence interval counterfactual exposure g(x) = 0.5x
#> $pif
#> [1] 0.1789588
#> 
#> $ci
#> [1] 0.06255514 0.29536242
```

Now suppose the user does not have access to the individual-level
exposure data but does know its mean and variance. This is typical in
publications, where only the mean and variance are reported. We can
approximate the PAF using only the first two moments.

``` r
paf.app(mean(x), var(x), length(x), log(1.27), 0.002)
```
