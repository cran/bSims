## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
par(mar = c(1, 1, 1, 1))
set.seed(429)
suppressPackageStartupMessages(library(bSims))

## ----sim-bin,eval=FALSE-------------------------------------------------------
#  p <- 0.5
#  
#  Y <- rbinom(1, size = 1, prob = p)

## ----sim-pois,eval=FALSE------------------------------------------------------
#  D <- 2 # individuals / unit area
#  A <- 1 # area
#  p <- 0.8 # probability of availability given presence
#  q <- 0.5 # probability of detection given availability
#  
#  N <- rpois(1, lambda = A * D)
#  Y <- rbinom(1, size = N, prob = p * q)

## ----intro,fig.width=15,fig.height=10,out.width='100%'------------------------
library(bSims)

phi <- 0.5                 # singing rate
tau <- 1:3                 # detection distances by strata
tbr <- c(3, 5, 10)         # time intervals
rbr <- c(0.5, 1, 1.5)      # count radii

l <- bsims_init(extent=10, # landscape
  road=0.25, edge=0.5)
p <- bsims_populate(l,     # population
  density=c(1, 1, 0))
e <- bsims_animate(p,      # events
  vocal_rate=phi,
  move_rate=1, movement=0.2)
d <- bsims_detect(e,       # detections
  tau=tau)
x <- bsims_transcribe(d,   # transcription
  tint=tbr, rint=rbr)

get_table(x) # removal table

op <- par(mfrow=c(2,3), cex.main=2)
plot(l, main="Initialize")
plot(p, main="Populate")
plot(e, main="Animate")
plot(d, main="Detect")
plot(x, main="Transcribe")
par(op)

## ----landscape,fig.width=10,fig.height=10,out.width='100%'--------------------
(l1 <- bsims_init(extent = 10, road = 0, edge = 0, offset = 0))
(l2 <- bsims_init(extent = 10, road = 1, edge = 0, offset = 0))
(l3 <- bsims_init(extent = 10, road = 0.5, edge = 1, offset = 2))
(l4 <- bsims_init(extent = 10, road = 0, edge = 5, offset = 5))

op <- par(mfrow = c(2, 2))
plot(l1, main = "Habitat")
points(0, 0, pch=3)
plot(l2, main = "Habitat & road")
lines(c(0, 0), c(-5, 5), lty=2)
plot(l3, main = "Habitat, edge, road + offset")
arrows(0, 0, 2, 0, 0.1, 20)
lines(c(2, 2), c(-5, 5), lty=2)
points(0, 0, pch=3)
plot(l4, main = "2 habitats")
arrows(0, 0, 5, 0, 0.1, 20)
lines(c(5, 5), c(-5, 5), lty=2)
points(0, 0, pch=3)
par(op)

## ----pop-pois-----------------------------------------------------------------
bsims_populate(l1)

## ----pop-nb-------------------------------------------------------------------
summary(rpois(100, 100)) # Poisson variation
summary(MASS::rnegbin(100, 100, 0.8)) # NegBin variation
negbin <- function(lambda, ...) MASS::rnegbin(1, lambda, ...)
bsims_populate(l1, abund_fun = negbin, theta = 0.8)
## constant abundance
bsims_populate(l1, abund_fun = function(lambda, ...) lambda)

## ----pop-xy,fig.height=9,fig.width=6------------------------------------------
D <- 0.5
f_abund <- function(lambda, ...) lambda

## systematic
f_syst <- function(d)
  (1-exp(-d^2/1^2) + dlnorm(d, 2)/dlnorm(exp(2-1),2)) / 2
## clustered
f_clust <- function(d)
  exp(-d^2/1^2) + 0.5*(1-exp(-d^2/4^2))

p1 <- bsims_populate(l1, density = D, abund_fun = f_abund)
p2 <- bsims_populate(l1, density = D, abund_fun = f_abund, xy_fun = f_syst)
p3 <- bsims_populate(l1, density = D, abund_fun = f_abund, xy_fun = f_clust)

distance <- seq(0,10,0.01)
op <- par(mfrow = c(3, 2))
plot(distance, rep(1, length(distance)), type="l", ylim = c(0, 1), 
  main = "random", ylab=expression(f(d)), col=2)
plot(p1)

plot(distance, f_syst(distance), type="l", ylim = c(0, 1), 
  main = "systematic", ylab=expression(f(d)), col=2)
plot(p2)

plot(distance, f_clust(distance), type="l", ylim = c(0, 1), 
  main = "clustered", ylab=expression(f(d)), col=2)
plot(p3)
par(op)

## ----pop-dens,fig.width=10,fig.height=10,out.width='100%'---------------------
D <- c(H = 2, E = 0.5, R = 0)

op <- par(mfrow = c(2, 2))
plot(bsims_populate(l1, density = D), main = "Habitat")
plot(bsims_populate(l2, density = D), main = "Habitat & road")
plot(bsims_populate(l3, density = D), main = "Habitat, edge, road + offset")
plot(bsims_populate(l4, density = D), main = "2 habitats")
par(op)

## ----beh-events,fig.height=4,fig.width=6--------------------------------------
l <- bsims_init()
p <- bsims_populate(l, density = 0.5)
e1 <- bsims_animate(p, vocal_rate = 1)

head(get_events(e1))
plot(get_events(e1))
curve((1-exp(-1*x)) * get_abundance(e1), col=2, add=TRUE)

## ----beh-move,fig.width=10,fig.height=5,out.width='100%'----------------------
e2 <- bsims_animate(p, move_rate = 1, movement = 0.25)

op <- par(mfrow = c(1, 2))
plot(e1, main = "Closure")
plot(e2, main = "Movement")
par(op)

## ----beh-mix,fig.width=6,fig.height=4-----------------------------------------
e3 <- bsims_animate(p, 
  vocal_rate = c(25, 1), mixture = c(0.33, 0.67))

plot(get_events(e3))
curve((1-0.67*exp(-1*x)) * get_abundance(e3), col=2, add=TRUE)

## ----beh-move2,fig.width=5,fig.height=5---------------------------------------
plot(bsims_animate(bsims_populate(l4, density = D), 
  move_rate = c(0.5, 1, 1), movement = c(0, 0.2, 0.2), 
  mixture = 1), main="Strata based mixtures")

## ----beh-move3,fig.width=10,fig.height=5,out.width='100%'---------------------
op <- par(mfrow = c(1, 2))
plot(bsims_animate(bsims_populate(l2, density = D), 
  move_rate = 1, movement = 0.3, 
  avoid = "none"), main="Movement not restricted")
plot(bsims_animate(bsims_populate(l2, density = D), 
  move_rate = 1, movement = 0.3, 
  avoid = "R"), main="Movement restricted")
par(op)

## ----beh-move4,fig.width=10,fig.height=5,out.width='100%'---------------------
e4 <- update(e2, allow_overlap=FALSE)

op <- par(mfrow = c(1, 2))
plot(e2, main = "Overlap")
plot(e2$tess, add=TRUE, wlines="tess",
  showpoints=FALSE, cmpnt_lty=1)
plot(e4, main = "No overlap")
plot(e4$tess, add=TRUE, wlines="tess",
  showpoints=FALSE, cmpnt_lty=1)
par(op)

## ----det,fig.width=5, fig.height=5--------------------------------------------
(d1 <- bsims_detect(e2, tau = 2))

head(get_detections(d1))
plot(d1)

## ----segm,fig.height=8,fig.width=6--------------------------------------------
tau <- c(1, 2, 3, 2, 1)
d <- seq(0, 4, 0.01)
dist_fun <- function(d, tau) exp(-(d/tau)^2) # half normal
#dist_fun <- function(d, tau) exp(-d/tau) # exponential
#dist_fun <- function(d, tau) 1-exp(-(d/tau)^-2) # hazard rate

b <- c(0.5, 1, 1.5, 2) #  boundaries

op <- par(mfrow=c(2, 1))
plot(d, dist_fun2(d, tau[1], dist_fun), type="n",
  ylab="g(d)", xlab="d (100 m)", axes=FALSE,
  main="Sound travels from left to right")
axis(1)
axis(2)
for (i in seq_len(length(b)+1)) {
  x1 <- c(0, b, 4)[i]
  x2 <- c(0, b, 4)[i+1]
  polygon(c(0, b, 4)[c(i, i, i+1, i+1)], c(0, 1, 1, 0),
    border=NA,
    col=c("darkolivegreen1", "burlywood1", "lightgrey",
    "burlywood1", "darkolivegreen1")[i])
}
lines(d, dist_fun2(d, tau[1], dist_fun))
lines(d, dist_fun2(d, tau[2], dist_fun))
lines(d, dist_fun2(d, tau[3], dist_fun))
lines(d, dist_fun2(d, tau, dist_fun, b), col=2, lwd=3)

plot(rev(d), dist_fun2(d, tau[1], dist_fun), type="n",
  ylab="g(d)", xlab="d (100 m)", axes=FALSE,
  main="Sound travels from right to left")
axis(1)
axis(2)
for (i in seq_len(length(b)+1)) {
  x1 <- c(0, b, 4)[i]
  x2 <- c(0, b, 4)[i+1]
  polygon(c(0, b, 4)[c(i, i, i+1, i+1)], c(0, 1, 1, 0),
    border=NA,
    col=c("darkolivegreen1", "burlywood1", "lightgrey",
    "burlywood1", "darkolivegreen1")[i])
}
lines(rev(d), dist_fun2(d, tau[1], dist_fun))
lines(rev(d), dist_fun2(d, tau[2], dist_fun))
lines(rev(d), dist_fun2(d, tau[3], dist_fun))
lines(rev(d), dist_fun2(d, tau, dist_fun, rev(4-b)), col=2, lwd=3)
par(op)

## ----det2,fig.width=5,fig.height=5--------------------------------------------
e5 <- bsims_animate(
  bsims_populate(
    bsims_init(road = 0.2, edge = 0.4), density = D), 
  move_rate = 1, movement = 0.2)
(d2 <- bsims_detect(e5, tau = c(1, 2, 3), event_type = "both"))

head(get_detections(d2))
plot(d2)

## ----det3,fig.width=5,fig.height=5--------------------------------------------
(d3 <- bsims_detect(e2, tau = c(1.5, 3), event_type = "both"))

dtab <- get_detections(d3)
tmp <- with(dtab, table(i, v))
c("heard"=sum(tmp[,"0"] == 0 & tmp[,"1"] > 0),
  "seen"=sum(tmp[,"0"] > 0 & tmp[,"1"] == 0),
  "heard & seen"=sum(tmp[,"0"] > 0 & tmp[,"1"] > 0))

plot(d3)

## ----trans,fig.height=6,fig.width=6-------------------------------------------
x <- bsims_transcribe(d1,
  tint = c(2, 4, 6, 8, 10), 
  rint = c(0.5, 1, 1.5, Inf))
x
plot(x)

get_table(x, "removal")
get_table(x, "visits")

## ----H,fig.height=5,fig.width=5-----------------------------------------------
phi <- 0.5 # singing rate
tau <- 2   # detection distance
Den <- 10  # density

set.seed(1)
l <- bsims_init()
a <- bsims_populate(l, density=Den)
b <- bsims_animate(a, vocal_rate=phi)
o <- bsims_detect(b, tau=tau)

tint <- c(1, 2, 3, 4, 5)
rint <- c(0.5, 1, 1.5, 2) # truncated at 200 m
(x <- bsims_transcribe(o, tint=tint, rint=rint))
(y <- get_table(x, "removal")) # binned new individuals
colSums(y)
rowSums(y)
plot(x)

## ----Dtr----------------------------------------------------------------------
library(detect)
cbind(true = c(phi=phi, tau=tau, D=Den), 
  estimate = estimate(x))

## ----Dinf---------------------------------------------------------------------
rint <- c(0.5, 1, 1.5, 2, Inf) # unlimited

(x <- bsims_transcribe(o, tint=tint, rint=rint))
(y <- get_table(x, "removal")) # binned new individuals
colSums(y)
rowSums(y)

cbind(true = c(phi=phi, tau=tau, D=Den), 
  estimate = estimate(x))

