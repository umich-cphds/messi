
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package `messi`

# Mediation Analysis with External Summary-Level Information on the Total Effect of Exposure

[![](https://img.shields.io/badge/devel%20version-0.1.1-blue.svg)](https://github.com/umich-cphds/messi)

## Overview

This `R` package fits the hard constraint, soft constraint, and
unconstrained models in Boss et al. (2023) for mediation analyses with
external summary-level information on the total effect.

## Installation

If the devtools package is not yet installed, install it first:

``` r
install.packages('devtools')
```

``` r
# install the package from Github:
devtools::install_github('umich-cphds/messi') 
```

Once installed, load the package:

``` r
library(messi)
```

## Example Usage

For this example, we simulate data and test the cases of null and
non-null mediation effect.

### Test case with null mediation effect

``` r
library(MASS)
set.seed(20230419)
n <- 500
p <- 20
Y <- rnorm(n = n, mean = 0, sd = 1)
M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
A <- rnorm(n = n, mean = 0, sd = 1)
C <- NULL
method <- "Unconstrained"
s2.fixed <- NULL
T.hat.external <- NULL
var.T.hat.external <- NULL

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

set.seed(20230419)
n <- 500
p <- 20
Y <- rnorm(n = n, mean = 0, sd = 1)
M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
A <- rnorm(n = n, mean = 0, sd = 1)
C <- NULL
method <- "Hard"
s2.fixed <- NULL
T.hat.external <- 0
var.T.hat.external <- NULL

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

set.seed(20230419)
n <- 500
p <- 20
Y <- rnorm(n = n, mean = 0, sd = 1)
M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
A <- rnorm(n = n, mean = 0, sd = 1)
C <- NULL
method <- "Soft EB"
s2.fixed <- NULL
T.hat.external <- 0
var.T.hat.external <- 0.2

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

set.seed(20230419)
n <- 500
p <- 20
Y <- rnorm(n = n, mean = 0, sd = 1)
M <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(1, p))
A <- rnorm(n = n, mean = 0, sd = 1)
C <- NULL
method <- "Soft Fixed"
s2.fixed <- 1
T.hat.external <- 0
var.T.hat.external <- 0.2

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
```

### Test case with non-null mediation effect

``` r
data(Med)

Y = Med$Y
M = Med$M
A = Med$A
C = Med$C
T.hat.external = Med$T.hat.external
var.T.hat.external = Med$var.T.hat.external

test <- messi(Y = Y, M = M, A = A, C = C, method = 'Unconstrained', T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = NULL)

n = Med$n
p = Med$p

plot_messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, 
           labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
              
test <- messi(Y = Y, M = M, A = A, C = C, method = 'Hard', T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = NULL)
```
