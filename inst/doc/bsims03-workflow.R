## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
par(mar = c(1, 1, 1, 1))
set.seed(429)
suppressPackageStartupMessages(library(bSims))

## ----sim1---------------------------------------------------------------------
library(bSims)

tint <- c(2, 4, 6, 8, 10)
rint <- c(0.5, 1, 1.5, 2, Inf) # unlimited

## no road
b1 <- bsims_all(
  road = 0,
  density = c(1, 1, 0),
  tint = tint,
  rint = rint)
## road
b2 <- bsims_all(
  road = 0.5,
  density = c(1, 1, 0),
  tint = tint,
  rint = rint)
b1
b2

## ----sim2---------------------------------------------------------------------
b1$new()
b2$new()

## ----sim3,eval=FALSE----------------------------------------------------------
#  B <- 25  # number of runs
#  bb1 <- b1$replicate(B)
#  bb2 <- b2$replicate(B)

## ----grid1--------------------------------------------------------------------
s <- expand_list(
  road = c(0, 0.5, 1),
  density = list(c(1, 1, 0)),
  tint = list(tint),
  rint = list(rint))
str(s)


## ----grid2,eval=FALSE---------------------------------------------------------
#  b <- lapply(s, bsims_all)
#  nc <- 4 # number of cores to use
#  library(parallel)
#  cl <- makeCluster(nc)
#  bb <- lapply(b, function(z) z$replicate(B, cl=cl))
#  stopCluster(cl)

## ----grid3--------------------------------------------------------------------
s <- expand_list(
  road = c(0, 0.5),
  xy_fun = list(
    NULL,
    function(d) exp(-d^2/1^2) + 0.5*(1-exp(-d^2/4^2))),
  density = list(c(1, 1, 0)),
  tint = list(tint),
  rint = list(rint))
str(s)

