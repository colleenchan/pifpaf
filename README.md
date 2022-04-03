
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
relative risk function takes the form
![\\text{RR}(x) = \\exp(x) = 1.27](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BRR%7D%28x%29%20%3D%20%5Cexp%28x%29%20%3D%201.27 "\text{RR}(x) = \exp(x) = 1.27")
with variance 0.002.

``` r
library(pifpaf)
# generate some exposure data
set.seed(1)
x <- rweibull(1000, 1.2, 1.66) 

# Estimate PAF
pif.ind(x, beta = log(1.27), varbeta = 0.002)
#> Estimating PAF and 95% confidence interval
#> $pif
#> [1] 0.3447482
#> 
#> $ci
#> [1] 0.2273862 0.4621102
```

For the PIF, we allow counterfactual exposures of the form
![g(x) = a + bx](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g%28x%29%20%3D%20a%20%2B%20bx "g(x) = a + bx").

``` r
# Estimate PIF for a counterfactual exposure of a 1 unit decrease
pif.ind(x, beta = log(1.27), varbeta = 0.002, a = -1)
#> Estimating PIF with counterfactual exposure g(x) = x - 1 and 95% confidence interval
#> $pif
#> [1] 0.198881
#> 
#> $ci
#> [1] 0.08201511 0.31574683

# Estimate PIF for a counterfactual exposure of 50% decrease
pif.ind(x, beta = log(1.27), varbeta = 0.002, b = 0.5)
#> Estimating PIF with counterfactual exposure g(x) = 0.5x and 95% confidence interval
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
paf.app(meanx = 1.55, varx = 1.6, n = 1000, beta = log(1.27), varbeta = 0.002)
```
