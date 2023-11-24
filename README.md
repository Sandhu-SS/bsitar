<!-- README.md is generated from README.Rmd. Please edit that file -->

## bsitar

The **bsitar** package provides an interface for Bayesian implementation
of Super Imposition by Translation and Rotation (SITAR) growth model.
The SITAR is a shape-invariant nonlinear mixed effect model that fits a
natural cubic spline mean curve and aligns individual-specific growth
curves to the underlying mean curve via a set of random effects: the
size, timing and intensity. The **bsitar** package package is a
front-end to the R package **brms** which uses the **Stan** program (see
<https://mc-stan.org/>) to performing full Bayesian inference for a
range of regression specifications including hierarchical multivariate
modeling.

## Installation

You can install the development version of bsitar from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Sandhu-SS/bsitar")
```

## Example

This is a basic example which shows you how to fit the model:

``` r
library(bsitar)
```

``` r
data(heights, package = "sitar")
data <- heights
rm(heights)

data <- data %>% 
  dplyr::select(id, age, height) %>% 
  tidyr::drop_na(height) %>% 
  dplyr::rename(y = height) %>% 
  dplyr::arrange(id, age) %>% 
  dplyr::mutate(id = as.factor(id) ) 
```

``` r
fit1 <- bsitar(age, y, id, data = data, df = 3, chains = 2, iter = 1000)
```
