
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package `messi`

# Mediation Analysis with External Summary-Level Information on the Total Effect of Exposure

[![](https://img.shields.io/badge/devel%20version-0.1.1-blue.svg)](https://github.com/umich-cphds/messi)

## Overview

This `R` package is an ….

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

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

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

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

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

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

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

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
```

### Test case with non-null mediation effect

``` r
#Generate Data from Mediation Model
set.seed(20220222)
n <- 200
p.con <- 5
p.mediators <- 50
rho.con.exp <- 0.2
rho.mediators <- 0.2
r2.mediator <- 0.05
total.effect.internal <- 1
r2.outcome <- 0.2
n.external.ratio <- 100
is.same.external <- TRUE
total.effect.external <- total.effect.internal
sim.dat <- sim.data(n = n, p.con = p.con, p.mediators = p.mediators, rho.con.exp = rho.con.exp,
                    rho.mediators = rho.mediators, r2.mediator = r2.mediator, r2.outcome = r2.outcome,
                    total.effect.internal = total.effect.internal,
                    n.external.ratio = n.external.ratio, is.same.external = is.same.external,
                    total.effect.external = total.effect.external)

#True parameter values that generated the data
alpha.a.true <- sim.dat$alpha.a
alpha.c.true <- sim.dat$alpha.c
Sigma.m.true <- sim.dat$Sigma.m
beta.m.true <- sim.dat$beta.m
beta.a.true <- sim.dat$beta.a
beta.c.true <- sim.dat$beta.c
sigma.e.sq.true <- sim.dat$sigma.e.sq

#Estimated Internal and External Total Effects
T.hat.internal <- sim.dat$T.hat.internal
var.T.hat.internal <- sim.dat$var.T.hat.internal
T.hat.external <- sim.dat$T.hat.external
var.T.hat.external <- sim.dat$var.T.hat.external

set.seed(20230419)
Y <- sim.dat$Y
M <- sim.dat$M
A <- sim.dat$A
C <- sim.dat$C
method <- "Unconstrained"
s2.fixed <- NULL

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

set.seed(20230419)
Y <- sim.dat$Y
M <- sim.dat$M
A <- sim.dat$A
C <- sim.dat$C
method <- "Hard"
s2.fixed <- NULL

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

set.seed(20230419)
Y <- sim.dat$Y
M <- sim.dat$M
A <- sim.dat$A
C <- sim.dat$C
method <- "Soft EB"
s2.fixed <- NULL

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)

set.seed(20230419)
Y <- sim.dat$Y
M <- sim.dat$M
A <- sim.dat$A
C <- sim.dat$C
method <- "Soft Fixed"
s2.fixed <- 1

test <- messi(Y = Y, M = M, A = A, C = C, method = method, T.hat.external = T.hat.external,
              var.T.hat.external = var.T.hat.external, s2.fixed = s2.fixed)

plot.messi(n = n, alpha.a.hat = test$alpha.a.hat, beta.m.hat = test$beta.m.hat, labels = paste0("M",1:p), asym.var.mat = test$asym.var.mat)
```