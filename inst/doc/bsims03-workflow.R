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

## ----equal--------------------------------------------------------------------
bsims_all(
  road = 0.5,
  density = 1)

bsims_all(
  list(
    road = 0.5,
    density = 1))

bsims_all(
  data.frame(
    road = 0.5,
    density = 1))

## ----variables----------------------------------------------------------------
# number of stations to visit
n <- 5

# random predictors: continuous and discrete
x <- data.frame(x1=runif(n,-1,2), x2=rnorm(n))

# density
D <- drop(exp(model.matrix(~x2, x) %*% c(0,-0.5)))
summary(D)

# cue rate
phi <- drop(exp(model.matrix(~x1+I(x1^2), x) %*% c(-1,-0.25,-1)))
summary(phi)

# this data frame collects the columns to be used as arguments
s <- data.frame(
    D=D,
    vocal_rate = phi, 
    duration = 10,
    condition = "det1",
    tau = 1)

# each row from s becomes a simulation settings object
bb <- lapply(1:n, function(i) bsims_all(s[i,]))

# define how you want the data extracted
get_counts <- function(b) {
    o <- b$new() # simulate
    get_table(o)[1,1]
}

x$y <- sapply(bb, get_counts)
x

