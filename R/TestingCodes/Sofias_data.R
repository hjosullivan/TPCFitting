library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(cores=3)

source('response_functions.R')

dat <- read.csv('Downloads/RichardAlexSubset.csv', sep=',')
# My functions expect 'value'.
setnames(dat, 'OriginalTraitValue', 'value')
dat$temp.K <- dat$ConTemp + 273.15
# Get rid of traits < 0.
dat <- dat[-which(dat$value <= 0), ]
dat.small <- dat[, c('OriginalID', 'temp.K', 'value')]
dat.small <- data.table(dat.small)

# Get rid of records with fewer than n values.
n <- 5
n.records <- dat.small[, .N, by=OriginalID]
ids.few   <- n.records$OriginalID[which(n.records$N >= n)]
dat.small <- dat.small[dat.small$Original %in% ids.few, ]

trefs <- -10:40

crazy.pars <- foreach (i = 1:length(trefs)) %dopar%
{
  sf.pars <- dat.small[, SchoolFit(dat=.SD, Tref.K=trefs[i] + 273.15), by=.(OriginalID)]
  sf.pars$Tref <- trefs[i]
  sf.pars <- list(sf.pars)
}


