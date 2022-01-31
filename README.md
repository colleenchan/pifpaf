
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pifpaf

<!-- badges: start -->
<!-- badges: end -->

The goal of pifpaf is to estimate potential impact fractions and
population attributable fractions using individual-level or aggregate
data.

## Installation

You can install the released version of pifpaf from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("pifpaf")
```

## Example

This is a basic example which shows you how to estimate the PAF when
individual-level exposure data is available to the user, where the
relative risk function takes the form RR(*x*) = exp (*x*) = 1.27.

``` r
library(pifpaf)
x <- rweibull(1000, 1.2, 1.66) # synthetic data
paf.emp(x, log(1.27), 0.002)
```

Now suppose the user does not have access to the individual-level
exposure data but does know its mean and variance. This is typical in
publications, where only the mean and variance are reported. We can
approximate the PAF using only the first two moments.

``` r
paf.app(mean(x), var(x), length(x), log(1.27), 0.002)
```