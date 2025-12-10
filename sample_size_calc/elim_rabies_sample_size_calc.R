rm(list = ls())
library(lme4)
library(GLMMmisc)  # https://github.com/pcdjohnson/GLMMmisc
library(parallel)

# power analysis for superiority trial of community-led delivery of thermotolerant
# rabies vaccine

# function to report messages (e.g. progress) from inside mclapply
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}

# function to adjust binomial GLMM for jensen's inequality 
jensen.adjust <-
  function(p, V, method = "mcculloch", inverse = FALSE) {
    stopifnot(!(method == "mcculloch" & inverse))
    Beta <- qlogis(p)
    if(method == "mcculloch") {
      return(plogis(Beta - 0.5 * V * tanh(Beta * (1 + 2 * exp(-0.5 * V))/6)))
    }
    if(method == "zeger") {
      if(inverse) 
        return(plogis(Beta * sqrt(256 * V / (75 * pi^2) + 1))) else
          return(plogis(Beta/sqrt(1 + ((16 * sqrt(3))/(15 * pi))^2 * V)))
    }
  }

# quick wald-z p-value from lmer fits
wald.p <- 
  function(fit, term) {
    t.val <- coef(summary(fit))[term, "t value"]
    2 * (1 - pnorm(abs(t.val)))
  }


# function to simulate variance variability
sim.var <-
  function(mu, v, min.out = 0, use.mu = TRUE) {
    if(use.mu) return(mu)
    out <- -1
    while(out < min.out) {
      s <-  v/mu
      a <- mu/s
      out <- rgamma(1, shape = a, scale = s)
    }
    out
  }

# load pilot data to supply number of dogs per household
census <- read.csv("data/postvac_census.csv", stringsAsFactors = FALSE) # HUMAN AND DOG Census
names(census) <- tolower(names(census))
head(census)
census$totdogs <- census$dogs + census$pups
census$vacc_totdogs <- census$vacc_dogs + census$vacc_pups
census <- droplevels(census[census$totdogs > 0.5, ])
census$unvacc = census$totdogs - census$vacc_totdogs
census$cov <- census$vacc_totdogs / census$totdogs

# which variance estimates to use? pre-trial (estimated from pilot data
# and used in the sample size calculation)
# or post-trial (estimated from the trial data, used to assess the 
# sample size calculation assumptions)?
prepost <- c("pretrial", "posttrial")[2]

# specify the assumptions
params <-
  cbind(
    pretrial = c(hh.var = 11.5, sv.var = 4.4, vill.var = 0.55, ward.var = 0.023,
                 dogs.per.hh = 2.5, hh.per.sv = 10),
    posttrial = c(hh.var = 3.81, sv.var = 2.36, vill.var = 0, ward.var = 0.53,
                  dogs.per.hh = 1.42, hh.per.sv = 9))

# design:
n.yr <- 1 # no of years the expt is run for
n.trt <- 2  # trt1 is gov-led, trt2 is com-led
n.ward.per.trt <- 56
n.vill.per.ward <- 1
n.subvill.per.vill <- 3
n.hh.per.subvill <- params["hh.per.sv", prepost] * 2 * 3 # to simulate 3 years (nb 20 hh results in 10 hh per year in the final data set)
# 18 (resulting in 9) equates to the number in the final data set

# vary the variances?
var.var <- FALSE

# variance estimates. trt1 is gov-led, trt2 is com-led
ward.plus.vill.var <- function(...) sim.var(mu = params["ward.var", prepost], v = 0.05, use.mu = !var.var)
# 0.023 is a point estimate from pilot data, variance of this estimate = 0.05

vill.var <- function(...) sim.var(mu = params["vill.var", prepost], v = 0.93, use.mu = !var.var)
# assume zero because trial estimate combines ward and village variances, because 
# only one village was surveyed per ward.
# 0.55 is a point estimate from pilot data, variance of this estimate = 0.93

subvill.var <- function(...) sim.var(mu = params["sv.var", prepost], v = 2.6, use.mu = !var.var)
# 4.4 based on pilot data (Maganga). variance of this estimate = 2.6.

hh.var <- function(...) sim.var(mu = params["hh.var", prepost], v = 1.3, use.mu = !var.var)
# 11.5 based on pilot data. variance of this estimate = 1.3

# year to year variance
yr.var <- 0

# using the census data gives 2.5 dogs on average
# add more from a poisson with mean lambda
# Update based on trial, which surveyed 1.42 dogs/HH
# use 0.42 because we will +1 to avoid HH with zero dogs
lambda <- 0

# or, increase the number of dogs with a multiplier
dog.mult <- params["dogs.per.hh", prepost]/2.5

# how many dogs will be observed
dog.mult * (2.5 + lambda) * n.trt * n.ward.per.trt * n.vill.per.ward * 
  n.subvill.per.vill * n.hh.per.subvill * n.yr

# create template data frame
dat <- 
  expand.grid(hh = 1:n.hh.per.subvill,
              subvill = 1:n.subvill.per.vill, 
              vill = 1:n.vill.per.ward, 
              ward = factor(1:n.ward.per.trt), 
              trt = 1:n.trt,
              yr = factor(1:n.yr))

dat$yr.cent <- (dat$yr == 2) - 0.5

table(dat$vill)
dat$ward <- factor(paste(dat$trt, dat$ward, sep = "."))
dat$vill <- factor(paste(dat$ward, dat$vill, sep = "."))
dat$subvill <- factor(paste(dat$vill, dat$subvill, sep = "."))
dat$hh <- factor(paste(dat$subvill, dat$hh, sep = "."))


# each year, each dog goes into one of 4 categories based on vacc status at 1 and 11 months:

#   vaccinated at 1 and 11
#   vaccinated at 1 not 11
#   vaccinated at 11 not 1
#   unvaccinated at 1 and 11

# coverage
cv <- 
  list(gov = c(v1v11 = NA, v1u11 = NA, u1v11 = NA, u1u11 = NA),
       com = c(v1v11 = NA, v1u11 = NA, u1v11 = NA, u1u11 = NA)) 

# team led 60% at 1 month -> ~40% at 11 months
# CL 60% at 1 month -> 40% at 11 months, then 50% of new dogs are vaccinated

# birth = death rate:
lifespan <- 26/12  # years
period <- 10/12 # years

# lifespan ~ exp(lambda = 1/lifespan)


# team-led arm

cov01.gov <- 0.6   # coverage at the start of the year
cov11.gov <- cov01.gov * exp(-(1/lifespan) * period) # coverage after a year of dog turnover

cv$gov["v1u11"] <- cov01.gov - cov11.gov
cv$gov["v1v11"] <- cov11.gov
cv$gov["u1v11"] <- 0
cv$gov["u1u11"] <- 1 - cv$gov["v1u11"] - cv$gov["v1v11"] - cv$gov["u1v11"]

sum(cv$gov)

# community arm

cov01.com <- 0.6   # coverage at the start of the year
cov11.com <- cov01.com * exp(-(1/lifespan) * period) # coverage after a year of dog turnover, without follow-up vax

# assume that 50% of new dogs vaccinated
cv$com["v1u11"] <- cov01.com - cov11.com - (cov01.com - cov11.com) / 2
cv$com["v1v11"] <- cov11.com + cv$com["v1u11"]  
# assume that 25% of dogs missed at 1 month can be vaccinated at a rate of 50%:
# i.e 25% are potentially accessible to follow-up vax, and 50% are reached by the vaccinator
cv$com["u1v11"] <- 0.25 * 0.5 * (1 - cov01.com)
cv$com["u1u11"] <- 1 - cv$com["v1u11"] - cv$com["v1v11"] - cv$com["u1v11"]
sum(cv$com)

# check category proportions sum to 1
sapply(cv, sum)

# in the gov-led arm, coverage at v1:
cv$gov["v1v11"] + cv$gov["v1u11"]
# in the gov-led arm, coverage at v2:
cv$gov["v1v11"] + cv$gov["u1v11"]

# in the com-led arm, coverage at v1:
cv$com["v1v11"] + cv$com["v1u11"]
# in the gov-led arm, coverage at v2:
cv$com["v1v11"] + cv$com["u1v11"]

# coverage at v1 and v11:
sapply(cv, function(x) sum(x[substr(names(x), 1, 2) == "v1"]))
sapply(cv, function(x) sum(x[substr(names(x), 3, 5) == "v11"]))


# optionally simulate null hypothesis for checking type I error
###cv$com <- cv$gov


# calculate mean coverage over the year in the two arms
meancov <- 
  apply(rbind(
    sapply(cv, function(x) sum(x[substr(names(x), 1, 2) == "v1"])),
    sapply(cv, function(x) sum(x[substr(names(x), 3, 5) == "v11"]))), 2, mean)

# what odds ratio does this imply
(meancov["com"] / (1 - meancov["com"])) / (meancov["gov"] / (1 - meancov["gov"]))


# how many simulations to run
nsim <- 1000#00

# significance level
alpha <- 0.05


# run simulations
start.time <- Sys.time()
simres.list <-
  mclapply(1:nsim, function(i) {
    
    # print progress
    if(i/nsim == round(i/nsim, 1)) message_parallel(paste0(100 * (i/nsim), "%"))
    
    # to start with, just simulate 1 year
    rand.V <- 
      c(ward = ward.plus.vill.var(), 
        vill = vill.var(), 
        subvill = subvill.var(), 
        hh = hh.var())
    
    if(n.yr > 1) rand.V <- c(yr = yr.var, rand.V)
    
    
    dat$n.orig <-
      round(dog.mult * sample(census$totdogs, nrow(dat), replace = TRUE) + rpois(nrow(dat), lambda = lambda))
    mean(dat$n.orig)
    sum(dat$n.orig)
    dat[, names(cv$gov)] <- NA
    head(dat)
    
    for(trt in 1:length(cv)) {
      # trt = 2
      dat$n <- dat$n.orig
      dat$v1v11[dat$trt == trt] <- 
        sim.glmm(design.data = dat[dat$trt == trt, ], 
                 fixed.eff = 
                   c("(Intercept)" = 
                       as.vector(qlogis(jensen.adjust(cv[[trt]]["v1v11"], V = sum(rand.V), method = "zeger", inverse = TRUE)))),
                 rand.V = rand.V,
                 distribution = "binomial")$response
      dat[sort(sample(nrow(dat), 20)), ]
      
      dat$n[dat$trt == trt] <- dat$n.orig[dat$trt == trt] - dat$v1v11[dat$trt == trt]
      dat[sort(sample(nrow(dat), 20)), ]
      
      dat$v1u11[dat$trt == trt] <- 
        sim.glmm(design.data = dat[dat$trt == trt, ], 
                 fixed.eff = c("(Intercept)" = 
                                 as.vector(qlogis(jensen.adjust(cv[[trt]]["v1u11"]/(1 - cv[[trt]]["v1v11"]), 
                                                                V = sum(rand.V), method = "zeger", inverse = TRUE)))),
                 rand.V = rand.V,
                 distribution = "binomial")$response
      sum(dat$v1u11[dat$trt == trt]) / sum(dat$n.orig[dat$trt == trt])
      
      dat$n[dat$trt == trt] <- dat$n.orig[dat$trt == trt] - dat$v1v11[dat$trt == trt] - dat$v1u11[dat$trt == trt]
      dat[sort(sample(nrow(dat), 20)), ]
      
      dat$u1u11[dat$trt == trt] <- 
        sim.glmm(design.data = dat[dat$trt == trt, ], 
                 fixed.eff = 
                   c("(Intercept)" = 
                       as.vector(qlogis(jensen.adjust(cv[[trt]]["u1u11"]/sum(cv[[trt]][c("u1v11", "u1u11")]), V = sum(rand.V), method = "zeger", inverse = TRUE)))),
                 rand.V = rand.V,
                 distribution = "binomial")$response
      dat[sort(sample(nrow(dat), 20)), ]
      
      dat$u1v11[dat$trt == trt] <- dat$n[dat$trt == trt] - dat$u1u11[dat$trt == trt]
      dat$n <- NULL
      dat[sort(sample(nrow(dat), 20)), ]
      
    }
    
    dat$cov01 <- (dat$v1v11 + dat$v1u11) / dat$n.orig
    dat$cov11 <- (dat$v1v11 + dat$u1v11) / dat$n.orig
    dat$covmean <- (dat$cov01 + dat$cov11) / 2
    
    dat$trtbin <- dat$trt - 1
    dat$trt <- NULL
    table(dat$trtbin)
    
    dim(dat)
    dat.long <- reshape(dat, varying = list(c("cov01", "cov11")), direction = "long")
    dim(dat.long)
    dat.long$cov <- dat.long$cov01
    dat.long$cov01 <- NULL
    dat.long$timebin <- dat.long$time - 1
    dat.long$time <- NULL
    
    
    dat.long[dat.long$hh == levels(dat.long$hh)[1], ]
    
    # to simulate sampling different HH at each visit, drop half
    # of the houses at 1 month and the other half at 11 months
    
    levels(dat.long$hh)
    
    dat.long$hh.no <-
      as.numeric(sapply(strsplit(as.character(dat.long$hh), "\\."), function(x) rev(x)[1]))
    dim(dat.long)
    dat.long <-
      dat.long[
        (dat.long$hh.no < (n.hh.per.subvill/2 + 0.5) & dat.long$timebin == 0) |
          (dat.long$hh.no > (n.hh.per.subvill/2 + 0.5) & dat.long$timebin == 1), ]
    dim(dat.long)
    
    table(dat.long$subvill, dat.long$timebin)
    head(dat.long, 20)
    
    dat.long[dat.long$hh == levels(dat.long$hh)[10], ]
    
    nlevels(dat.long$hh) == length(dat.long$hh)
    
    dat.long$timebin.cent <- dat.long$timebin - 0.5
    
    
    form <-
      c(formmean = "covmean ~ trtbin + (1 | ward) + (1 | subvill)",
        form01 = "cov01 ~ trtbin + (1 | ward) + (1 | subvill)",
        form11 = "cov11 ~ trtbin + (1 | ward) + (1 | subvill)",
        formmean.diff = "cov ~ trtbin*timebin.cent + (1 | ward) + (1 | subvill)")
    
    if(n.yr > 1) form[1:4] <- paste(form[1:4], "+ yr.cent")
    
    # for speed, use the normal approximation to the binomial
    fitmean <- 
      lmer(form["formmean"], data = dat.long, weights = n.orig)
    fitmean.diff <- 
      lmer(form["formmean.diff"], data = dat.long, weights = n.orig)
    
    c(pwrmean.diff = wald.p(fitmean.diff, "trtbin") < alpha,
      estmean = cumsum(as.vector(fixef(fitmean)[c("(Intercept)", "trtbin")])),
      ci.width = as.vector(apply(confint(fitmean.diff, "trtbin", method = "Wald"), 1, diff)))
    
  }, mc.cores = detectCores())

finish.time <- Sys.time()
print(finish.time - start.time)

(res <- apply(do.call("rbind", simres.list), 2, mean))
