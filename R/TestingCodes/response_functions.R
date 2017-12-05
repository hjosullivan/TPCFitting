###############################################################################
# Temperature response functions.
# ~~~
# 2015-10-21
#
# Library of response functions, together with nonlinear fitting, to be used in
# data analysis code.
###############################################################################

BOLTZ <- 8.617e-5  # Boltzmann constant, in eV/K.



.GetE <- function(tmp.K, rate, Tp.K, k=BOLTZ)
{
  # Estimate starting value for E, taking linear regression using the rise part
  # of the curve only.
  # ~~~ Parameters ~~~
  # tmp.K : temperature data (in K).
  # rate  : rate data corresponding to temperature above.
  # Tp.K  : temperature at which rate peaks, used as a cutoff point.
  # k     : Boltzmann constant.

  tmp.K <- as.numeric(tmp.K)
  rate  <- as.numeric(rate)

  tmp.w <- which(tmp.K <= Tp.K)
  if (length(tmp.w) > 1)
  {
    m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp.K[tmp.w]))))
    return(abs(summary(m)$coefficients[2, c('Estimate', 'Std. Error')]))
  } else
  {
    return(c(0.7, 2))  # Arbitrary estimate if we can't do regression.
  }
}



.GetLnB0 <- function(tmp.K, rate, Tref.K)
{
  # Estimate starting value for the normalising constant.
  # ~~~ Parameters ~~~
  # tmp.K  : temperature data (in K).
  # rate   : rate data corresponding to temperature above.
  # Tref.K : estimate normalising constant at this temperature (in K).

  tmp.K <- as.numeric(tmp.K)
  rate  <- as.numeric(rate)

  if (min(tmp.K, na.rm=TRUE) > Tref.K)
  {
    return(log(min(rate[1], na.rm=TRUE)))
  } else
  {
    return(log(max(rate[which(tmp.K <= Tref.K)], na.rm=TRUE)))
  }
}



.GetTpk <- function(tmp.K, rate)
{
  # Temperature at which the rate is maximised (estimate of T.peak).
  # ~~~ Parameters ~~~
  # tmp.K : Temperature data (in K).
  # rate  : Rate data corresponding to temperature above.

  tmp.K <- as.numeric(tmp.K)
  rate  <- as.numeric(rate)

  return(max(tmp.K[which.max(rate)]))
}



Arr <- function(lnB0, E, Tref.K, tmp.K, k=BOLTZ)
{
  # (Boltzmann)-Arrhenius response function, in ln scale. 
  # ~~~ Parameters ~~~
  # lnB0   : Ln rate at T.ref.
  # E      : Activation energy.
  # Tref.K : Reference temperature at which the rate is 'lnB0' (in K).

  return(lnB0 - E/k * (1/tmp.K - 1/Tref.K))
}



ArrFit <- function(dat, trait.name, temp.name, Tref.K, k=BOLTZ)
{
  # Fit the Arrhenius model to the rising portion of the curve. This is mostly
  # for corroboration of the Schoolfield-Sharpe model, and also for estimates
  # for curves for which fitting failed.
  # ~~~ Parameters ~~~
  # dat            : data frame - with expected hardcoded column names.
  # trait.col.name : string -column name in 'dat' that corresponds to the trait.
  # temp.col.name  : string - column name in 'dat' that corresponds to
  #                  temperature (in K!).
  # Tref.K         : numeric - reference temperature (in K).
  # k              : numeric - Boltzmann constant.
  library(minpack.lm)
  library(sme)

  trait  <- dat[[trait.name]]
  temp.K <- dat[[temp.name]]

  if (min(temp.K, na.rm=TRUE) < 200)
  {
    stop('Error in Arrfit : make sure temperatures are in K')
  }

  # Using T.peak as cut off point, which means that the estimates will be a
  # little off. However, cutting off at a lower temperature than T.pk leads to
  # estimates of E that are substantially higher than those of the
  # Schoolfield-Sharpe model.
  Tpk.K  <- .GetTpk(tmp.K=temp.K, rate=trait)
  E.st   <- as.numeric(.GetE(tmp.K=temp.K, rate=trait, Tp.K=Tpk.K)[1])
  lnB.st <- .GetLnB0(tmp.K=temp.K, rate=trait, Tref.K=Tpk.K)
  dat    <- dat[temp.K < Tpk.K, ]

  arr.fit <- try(nlsLM(log(trait) ~ Arr(lnB0, E, Tref.K=Tref.K, tmp.K=temp.K),
                       start=c(lnB0=lnB.st, E=E.st),
                       lower=c(lnB0=-Inf, E=0),
                       upper=c(lnB0=Inf,  E=Inf),
                       control=list(minFactor=1 / 2^16, maxiter=1024)),
                 silent=TRUE)

  if (class(arr.fit) == 'try-error')
  {
    return(as.list(as.numeric(c(lnB0.st=NA, E.st=NA, lnB0=NA, E=NA, lnB0.se=NA,
                                E.se=NA, aic=NA, aic.c=NA))))
  } else
  {
    out <- as.list(summary(arr.fit)$coefficients[, c('Estimate', 'Std. Error')])
    out <- c(lnB.st, E.st, out, AIC(arr.fit), AICc(arr.fit))
    names(out) <- c('lnB0st', 'E.st', paste0(rep(names(coefficients(arr.fit)), 2), 
                                             rep(c('', '.se'), each=2)), 'aic', 'aic.c')
    return(out)
  }
}



School <- function(lnB0, E, E.D, Tp.K, Tref.K, tmp.K, k=BOLTZ, Tpk.form=FALSE)
{
  # Modified Schoolfield-Sharpe model (Schoolfield & Sharpe, 1981), used by Dan
  # and Gab in their E Lett paper.
  # ~~~ Parameters ~~~
  # lnB0     : Rate at T.ref.
  # E        : Activation energy.
  # E.D      : De-activation energy.
  # Tp.K     : T at which enzyme is 50% active and 50% high-temperature suppressed
  #            (in K), if Tpk.form=FALSE, else, temperature for max trait.
  # Tref.K   : Reference temperature at which the rate is ln.const (in K).
  # k        : Boltzmann constant.
  # tmp.K    : temperatures (in K), as value or vector.
  # Tpk.form : boolean - Fit Tpeak (Dimitris and Samraat) version (TRUE) or Th
  #            version (FALSE)?

  if (!Tpk.form)
  {
    return(lnB0 - E/k * (1/tmp.K - 1/Tref.K) - 
           log(1 + exp((E.D/k) * (1/Tp.K - 1/tmp.K))))
  } else
  {
    return(lnB0 - E/k * (1/tmp.K - 1/Tref.K) - 
           log(1 + (E / (E.D-E)) * exp((E.D/k * (1/Tp.K - 1/tmp.K)))))
  }
}



SchoolFit <- function(dat, trait.name, temp.name, Tref.K, Tpk.form=FALSE,
                      rand.st=FALSE, n.rand=100, k=BOLTZ)
{
  # Fit the Schoolfield model, and return NAs if no convergence. Optionally, use
  # randomised starting values, to improve chances of finding a fit.
  # ~~~ Parameters ~~~
  # dat        : data.frame or data.table - including temperature and trait data.
  # trait.name : string - name of column that holds trait values.
  # temp.name  : string - name of column that holds temperature data (in K).
  # Tref.K     : reference temperature (in K).
  # rand.st    : boolean - use random starting values?
  # n.rand     : numeric - number of random starting values used (if
  #              `rand.st=TRUE`).
  # k          : Boltzmann constant.
  library(minpack.lm)
  library(sme)
  library(truncnorm)

  trait  <- dat[[trait.name]]
  temp.K <- dat[[temp.name]]

  if (min(temp.K, na.rm=TRUE) < 200)
  {
    stop('Error in SchoolFit : make sure temperatures are in K')
  }

  # Estimate T.h as being approximately T.peak.
  Th.st  <- .GetTpk(tmp.K=temp.K, rate=trait)
  # This is now two values, one for slope, the other for standard error around
  # slope.
  E.st   <- as.numeric(.GetE(tmp.K=temp.K, rate=trait, Tp.K=Th.st))
  lnB.st <- .GetLnB0(tmp.K=temp.K, rate=trait, Tref.K=Tref.K)

  if (rand.st)
  {
    # Create randomised starting points.
    E.st.pe <- E.st[1]  # Slope value.
    Th.st  <- c(Th.st, rnorm(n.rand-1, mean=Th.st, sd=15))
    # We need truncated normal to ensure we don't get negative values of E.
    E.st   <- c(E.st.pe, rtruncnorm(n.rand-1, a=0, b=Inf, mean=E.st[1], sd=2 * E.st[2]))
    B.st   <- exp(lnB.st)
    # Randomise on linear scale. Again, we don't want negative rates.
    lnB.st <- c(lnB.st, log(rtruncnorm(n.rand-1, a=0, b=Inf, mean=B.st, sd=B.st / 2)))

    # We'll select the best model using AICc. Many of these turn out to be
    # similar.
    aics.out <- rep(NA, n.rand)

    for (i in 1:n.rand)
    {
      school.fit <- try(nlsLM(log(trait) ~ School(lnB0, E, E.D, Tp.K, Tref.K=Tref.K,
                                                  tmp.K=temp.K, Tpk.form=Tpk.form),
                              start=c(lnB0=lnB.st[i], E=E.st[i], E.D=4*E.st[i],
                                      Tp.K=Th.st[i]),
                              lower=c(lnB0=-Inf, E=0,  E.D=0,  Tp.K=250),
                              upper=c(lnB0=Inf,  E=30, E.D=50, Tp.K=350),
                              control=list(minFactor=1 / 2^16, maxiter=1024)),
                        silent=TRUE)

      if (class(school.fit) != 'try-error')
      {
        aics.out[i] <- AICc(school.fit)
      }
    }

    w <- which.min(aics.out)

    if (length(w) > 0)
    {
      school.fit <- try(nlsLM(log(trait) ~ School(lnB0, E, E.D, Tp.K, Tref.K=Tref.K,
                                                  tmp.K=temp.K, Tpk.form=Tpk.form),
                              start=c(lnB0=lnB.st[w], E=E.st[w], E.D=4*E.st[w],
                                      Tp.K=Th.st[w]),
                              lower=c(lnB0=-Inf, E=0,  E.D=0,  Tp.K=250),
                              upper=c(lnB0=Inf,  E=30, E.D=50, Tp.K=350),
                              control=list(minFactor=1 / 2^16, maxiter=1024)),
                        silent=TRUE)

      out <- as.list(summary(school.fit)$coefficients[, c('Estimate', 'Std. Error')])
      out <- c(lnB.st[w], E.st[w], 4*E.st[w], Th.st[w], out, AIC(school.fit),
               AICc(school.fit))
      names(out) <- c(paste0(rep(names(coefficients(school.fit)), 2), 
                             rep(c('.st', '', '.se'), each=4)), 'aic', 'aic.c')
      return(out)
    } else
    {
      return(as.list(as.numeric(c(lnB0.st=NA, E.st=NA, E.D.st=NA, Tp.K.st=NA,
                                  lnB0=NA,    E=NA,    E.D=NA,    Tp.K=NA,
                                  lnB0.se=NA, E.se=NA, E.D.se=NA, Tp.K.se=NA,
                                  aic=NA, aic.c=NA))))
    }
  } else
  {
    E.st <- as.numeric(E.st[1])  # Slope.
    school.fit <- try(nlsLM(log(trait) ~ School(lnB0, E, E.D, Tp.K, Tref.K=Tref.K,
                                                tmp.K=temp.K, Tpk.form=Tpk.form),
                            start=c(lnB0=lnB.st, E=E.st, E.D=4*E.st, Tp.K=Th.st),
                            lower=c(lnB0=-Inf,   E=0,    E.D=0,      Tp.K=0),
                            upper=c(lnB0=Inf,    E=Inf,  E.D=Inf,    Tp.K=Inf),
                            control=list(minFactor=1 / 2^16, maxiter=1024)),
                      silent=TRUE)

    if (class(school.fit) == 'try-error')
    {
      return(as.list(as.numeric(c(lnB0.st=NA, E.st=NA, E.D.st=NA, Tp.K.st=NA,
                                  lnB0=NA,    E=NA,    E.D=NA,    Tp.K=NA,
                                  lnB0.se=NA, E.se=NA, E.D.se=NA, Tp.K.se=NA,
                                  aic=NA, aic.c=NA))))
    } else
    {
      out <- as.list(summary(school.fit)$coefficients[, c('Estimate', 'Std. Error')])
      out <- c(lnB.st, E.st, 4*E.st, Th.st, out, AIC(school.fit), AICc(school.fit))
      names(out) <- c(paste0(rep(names(coefficients(school.fit)), 2), 
                             rep(c('.st', '', '.se'), each=4)), 'aic', 'aic.c')
      return(out)

      # We'll want something like this eventually, with variance-covariance
      # matrix.
      #return(cbind(summary(school.fit)$coefficients[, c('Estimate', 'Std. Error')],
      #vcov(school.fit)))
    }
  }
}



SchoolFitRes <- function(dat, trait.name, temp.name, Tref.K, Tpk.form=FALSE,
                         k=BOLTZ, res=FALSE)
{
  # Fit the Schoolfield model, return residuslas, and return NAs if no convergence.
  # ~~~ Parameters ~~~
  # dat        : data.frame or data.table - including temperature and trait data.
  # trait.name : string - name of column that holds trait values.
  # temp.name  : string - name of column that holds temperature data (in K).
  # Tref.K     : reference temperature (in K).
  # k          : Boltzmann constant.
  # res        : Boolean - return residuals instead of coefficients?
  library(minpack.lm)
  library(sme)
  library(nlstools)

  trait  <- dat[[trait.name]]
  temp.K <- dat[[temp.name]]

  if (min(temp.K, na.rm=TRUE) < 200)
  {
    stop('Error in SchoolFit : make sure temperatures are in K')
  }

  # Estimate T.h as being approximately T.peak.
  Th.st  <- .GetTpk(tmp.K=temp.K, rate=trait)
  E.st   <- .GetE(tmp.K=temp.K, rate=trait, Tp.K=Th.st)
  lnB.st <- .GetLnB0(tmp.K=temp.K, rate=trait, Tref.K=Tref.K)

  school.fit <- try(nlsLM(log(trait) ~ School(lnB0, E, E.D, Tp.K, Tref.K=Tref.K,
                                              tmp.K=temp.K, Tpk.form=Tpk.form),
                          start=c(lnB0=lnB.st, E=E.st, E.D=4*E.st, Tp.K=Th.st),
                          lower=c(lnB0=-Inf,   E=0,    E.D=0,      Tp.K=0),
                          upper=c(lnB0=Inf,    E=Inf,  E.D=Inf,    Tp.K=Inf),
                          control=list(minFactor=1 / 2^16, maxiter=1024)),
                    silent=TRUE)

  if (class(school.fit) == 'try-error')
  {
    return(list(`Fitted values`=NA, `Standardized residuals`=NA))
  } else
  {
    res <- nlsResiduals(school.fit)
    out <- list(`Fitted values`=res[[2]][, 1],
                `Standardized residuals`=res[[2]][, 2])
  }
}



SchoolLowIn <- function(lnB0, E, E.l.D, E.h.D, Th.l.K, Th.h.K, Tref.K, tmp.K,
			k=BOLTZ)
{
  # Modified Schoolfield-Sharpe model (Schoolfield & Sharpe, 1981), including
  # low temperature inactivation.
  # ~~~ Parameters ~~~
  # lnB0   : Rate at T.ref.
  # E      : Activation energy.
  # E.h.D  : High temperature de-activation energy.
  # E.l.D  : Low temperature de-activation energy.
  # Th.h.K : Low temperature Th.K.
  # Th.l.K : Low temperature Th.K.
  # Tref.K : Reference temperature at which the rate is ln.const (in K).
  # tmp.K  : temperatures (in K), as value or vector.
  # k      : Boltzmann constant.

  return(lnB0 - E/k * (1/tmp.K - 1/Tref.K) - 
         log(1 + exp((E.h.D/k) * (1/Th.h.K - 1/tmp.K)) + 
             exp((E.l.D/k) * (1/Th.l.K - 1/tmp.K))))
}



SchoolLowInFit <- function(dat, trait.name, temp.name, Tref.K, k=BOLTZ)
{
  # Fit the Schoolfield model, and return NAs if no convergence.
  # ~~~ Parameters ~~~
  # dat        : data.frame or data.table - including temperature and trait data.
  # trait.name : string - name of column that holds trait values.
  # temp.name  : string - name of column that holds temperature data (in K).
  # Tref.K     : reference temperature (in K).
  # k          : Boltzmann constant.
  library(minpack.lm)
  library(sme)

  trait  <- dat[[trait.name]]
  temp.K <- dat[[temp.name]]

  if (min(temp.K, na.rm=TRUE) < 200)
  {
    stop('Error in SchoolLowInFit : make sure temperatures are in K')
  }

  # Estimate T.h as being approximately T.peak.
  Th.st  <- .GetTpk(tmp.K=temp.K, rate=trait)
  E.st   <- .GetE(tmp.K=temp.K, rate=trait, Tp.K=Th.st)
  lnB.st <- .GetLnB0(tmp.K=temp.K, rate=trait, Tref.K=Tref.K)

  # This is a shit model and only fits like shit with port algorithm, or else it
  # doesn't fit at all.
  s.fit <- try(nlsLM(log(trait) ~ SchoolLowIn(lnB0, E, E.l.D, E.h.D, Th.l.K,
                                              Th.h.K, Tref.K=Tref.K,
                                              tmp.K=temp.K),
                     start=c(lnB0=lnB.st, E=E.st, E.h.D=4*E.st, E.l.D=-3*E.st,
                             Th.l.K=Th.st-10, Th.h.K=Th.st),
                     lower=c(lnB0=-Inf, E=0,   E.h.D=0,   E.l.D=-Inf, Th.l.K=0,   Th.h.K=0),
                     upper=c(lnB0=Inf,  E=Inf, E.h.D=Inf, E.l.D=0,    Th.l.K=Inf, Th.h.K=Inf),
                     data=dat, control=list(minFactor=1 / 2^16, maxiter=1024)),
               silent=TRUE)

  if (class(s.fit) == 'try-error')
  {
    return(as.list(as.numeric(c(lnB0.st=NA, E.st=NA, E.h.D.st=NA, E.l.D.st=NA,
                                Th.l.K.st=NA, Th.h.K.st=NA, 
                                lnB0=NA, E=NA, E.h.D=NA, E.l.D=NA, Th.l.K=NA,
                                Th.h.K=NA, lnB0.se=NA, E.se=NA, E.h.D.se=NA,
                                E.l.D.se=NA, Th.l.K.se=NA, Th.h.K.se=NA, aic=NA,
                                aic.c=NA))))
  } else
  {
    out <- as.list(summary(s.fit)$coefficients[, c('Estimate', 'Std. Error')])
    out <- c(lnB.st, E.st, 4*E.st, -3*E.st, Th.st-10, Th.st, 
             out, AIC(s.fit), AICc(s.fit))
    names(out) <- c(paste0(rep(names(coefficients(s.fit)), 2), 
                           rep(c('.st', '', '.se'), each=6)), 'aic', 'aic.c')
    return(out)
  }
}



Poly <- function(a1, a2, tpeak, alpha, p, temp)
{
  # Model from my paper with Dan, where
  # y = a1 * (-x + b)^alpha for x <= b, and
  # y = a2 * (x - b)^alpha for x >= b.
  # ~~~ Parameters ~~~
  # a1    : numeric - control rate of increase up to b.
  # a2    : numeric - control rate of decline after b.
  # tpeak : numeric - temperature at which rate peaks.
  # alpha : numeric - controls shape of the curve.
  # p     : numeric - peak rate.
  # temp  : numeric vector - vector of temperatures (units consistent with
  #         tpeak).
  
  x1 <- temp[temp <= tpeak]
  x2 <- temp[temp > tpeak]
  return(c(a1 * (-x1 + tpeak)^alpha + p, a2 * (x2 - tpeak)^alpha + p))
}



PolyFit <- function(dat, trait.name, temp.name)
{
  # Fit arbitrary pseudo-polynomial model (taken from my paper with Dan) to the
  # data, as a means of seeing if a completely alternative model can compete
  # with the Arrhenius/Schoolfield above.
  # ~~~
  # NOTE: Although this does ok in some cases, it often fits like shit. Yet, the
  # AICs are comparable to models above. Maybe check AIC_c? In any case, I think
  # this might well be pointless goose chase. It does not improve on
  # Schoolfield, it does not address Arrhenius breakpoints, it does not have
  # fewer parameters...
  library(minpack.lm)
  library(sme)

  trait  <- dat[[trait.name]]
  temp.K <- dat[[temp.name]]

  if (min(temp.K, na.rm=TRUE) < 200)
  {
    stop('Error in PolyFit : make sure temperatures are in K')
  }

  tpeak <- .GetTpk(tmp.K=temp.K, rate=trait)
  bmax  <- max(log(trait))

  poly.fit <- try(nlsLM(log(trait) ~ Poly(a1, a2, tpeak, alpha, p, temp=temp.K),
                        start=c(a1=-0.012, a2=-0.15, alpha=1.5, tpeak=tpeak, p=bmax),
                        lower=c(a1=-Inf,   a2=-Inf,  alpha=0.5, b=100,       p=-Inf),
                        upper=c(a1=0,      a2=0,     alpha=4,   b=400,       p=Inf),
                        data=dat, control=list(minFactor=1 / 2^16, maxiter=1e4)),
                  silent=TRUE)

  if (class(poly.fit) == 'try-error')
  {
    return(as.list(as.numeric(c(a1=NA, a2=NA, b=NA, alpha=NA, p=NA, aic=NA,
                                aic.c=NA))))
  } else
  {
    out <- as.list(summary(poly.fit)$coefficients[, 'Estimate'])
    out <- c(out, AIC(poly.fit), AICc(poly.fit))
    names(out) <- c(names(coefficients(poly.fit)), 'aic', 'aic.c')
    return(out)
  }
}
