###############################################################################
# Temperature response functions.
# ~~~
# 2015-10-21
#
# Library of response functions, together with nonlinear fitting, to be used in
# data analysis code.
###############################################################################

BOLTZ = 8.617e-5  # Boltzmann constant, in eV/K.



.GetE = function(tmp_K, rate, Tp_K, k=BOLTZ)
{
  # Estimate starting value for E, taking linear regression using the rise part
  # of the curve only.
  # ~~~ Parameters ~~~
  # tmp_K :numeric: temperature data (in K).
  # rate  :numeric: rate data corresponding to temperature above.
  # Tp_K  :numeric: temperature at which rate peaks, used as a cutoff point.
  # k     :numeric: Boltzmann constant.

  tmp_K = as.numeric(tmp_K)
  rate  = as.numeric(rate)

  tmp_w = which(tmp_K <= Tp_K)
  if (length(tmp_w) > 1)
  {
    m = lm(log(rate[tmp_w]) ~ I(1 / (k * (tmp_K[tmp_w]))))
    return(abs(summary(m)$coefficients[2, c('Estimate', 'Std. Error')]))
  } else
  {
    return(c(0.7, 2))  # Arbitrary estimate if we can't do regression.
  }
}



.GetLnB0 = function(tmp_K, rate, Tref_K)
{
  # Estimate starting value for the normalisation constant B0.
  # ~~~ Parameters ~~~
  # tmp_K  :numeric: temperature data (in K).
  # rate   :numeric: rate data corresponding to temperature above.
  # Tref_K :numeric: estimate normalising constant at this temperature (in K).

  tmp_K = as.numeric(tmp_K)
  rate  = as.numeric(rate)

  if (min(tmp_K, na.rm=TRUE) > Tref_K)
  {
    return(log(rate[1]))
  } else
  {
    return(log(max(rate[which(tmp_K <= Tref_K)], na.rm=TRUE)))
  }
}



.GetTpk = function(tmp_K, rate)
{
  # Temperature at which the rate is maximised (estimate of T.peak).
  # ~~~ Parameters ~~~
  # tmp_K :numeric: Temperature data (in K).
  # rate  :numeric: Rate data corresponding to temperature above.

  tmp_K = as.numeric(tmp_K)
  rate  = as.numeric(rate)

  return(max(tmp_K[which.max(rate)]))
}



Arr = function(lnB0, E, Tref_K, tmp_K, k=BOLTZ)
{
  # (Boltzmann)-Arrhenius response function, in ln scale. 
  # ~~~ Parameters ~~~
  # lnB0   :numeric: Ln rate at T_ref.
  # E      :numeric: Activation energy (eV).
  # Tref_K :numeric: Reference temperature at which the rate is 'lnB0' (K).

  return(lnB0 - E/k * (1/tmp_K - 1/Tref_K))
}



ArrFit = function(dat, trait_name, temp_name, unc_name=NULL, Tref_K, n_rand=1,
                  k=BOLTZ)
{
  # Fit the Arrhenius model to the rising portion of the curve. This is mostly
  # for corroboration of the Schoolfield-Sharpe model, and also for estimates
  # for curves for which fitting failed.
  # ~~~ Parameters ~~~
  # dat        :data.frame: with expected hardcoded column names.
  # trait_name :string: column name in 'dat' that corresponds to the trait.
  # temp_name  :string: column name in 'dat' that corresponds to
  #             temperature (in K).
  # unc_name   :string: column name in 'dat' with uncertainty to use to weight
  #             values in the fitting process.
  # Tref_K     :numeric: reference temperature (in K).
  # n_rand     :numeric: number of random starting values to use.
  # k          :numeric: Boltzmann constant.
  # ~~~ Output ~~~
  # Matrix with:
  # - Point estimates.
  # - Standard errors.
  # - Starting values used in the best fit.
  # - Variance covarience (2x2) matrix.
  # - AICc, AIC, BIC, and log likelihood values for teh fit.

  colnames(out) = c('pars', 'est', 'se', 'start', 'lnB0', 'E', 'AIC_c', 'AIC',
                    'BIC', 'log_Lik')
  library(minpack.lm)
  library(sme)
  library('truncnorm')


  trait  = dat[[trait_name]]
  temp_K = dat[[temp_name]]

  if (min(temp_K, na.rm=TRUE) < 200)
  {
    stop('Error in Arrfit : make sure temperatures are in K')
  }

  if (length(unc_name) > 0)
  {
    wgts = 1 / dat[[unc_name]]
  } else
  {
    l    = length(trait)
    wgts = rep(1 / l, l)
  }

  # This is two values: the slope, and the standard error around slope.
  E_st   = as.numeric(.GetE(tmp_K=temp_K, rate=trait, Tp_K=temp_K[which.max(trait)]))
  lnB_st = .GetLnB0(tmp_K=temp_K, rate=trait, Tref_K=Tref_K)
  B_st   = exp(lnB_st)

  # Create randomised starting points.
  E_st_pe = E_st[1]  # Slope value point estimate.

  if (n_rand > 1)
  {
    # We need truncated normal to ensure we don't get negative values of E.
    E_st = c(E_st_pe, rtruncnorm(n_rand-1, a=0, b=Inf, mean=E_st[1], sd=2 * E_st[2]))
    # Randomise on linear scale. Again, we don't want negative rates.
    lnB_st = c(lnB_st, log(rtruncnorm(n_rand-1, a=0, b=Inf, mean=B_st, sd=B_st / 2)))
  } else
  {
    E_st = E_st_pe
  }

  # We'll select the best model using AICc. Many of these are similar.
  aics_out = rep(NA, n_rand)

  for (i in 1:n_rand)
  {
    arr_fit = try(nlsLM(log(trait) ~ Arr(lnB0, E, Tref_K=Tref_K, tmp_K=temp_K),
                        start=c(lnB0=lnB_st[i], E=E_st[i]),
                        lower=c(lnB0=-Inf, E=0),
                        upper=c(lnB0=Inf,  E=30),
                        weights=wgts,
                        control=list(minFactor=1 / 2^16, maxiter=1024)),
                  silent=TRUE)

    if (class(arr_fit) != 'try-error')
    {
      aics_out[i] = AICc(arr_fit)
    }
  }

  w = which.min(aics_out)

  if (length(w) > 0 & n_rand == 1)
  {
    out = cbind(summary(arr_fit)$coefficients[, c('Estimate', 'Std. Error')],
                start=c(lnB_st, E_st), vcov(arr_fit), AIC_c=aics_out,
                AIC=AIC(arr_fit), BIC=BIC(arr_fit), log_Lik=logLik(arr_fit))
  } else if (length(w) > 0 & n_rand > 1)
  {
    arr_fit = try(nlsLM(log(trait) ~ Arr(lnB0, E, Tref_K=Tref_K, tmp_K=temp_K),
                        start=c(lnB0=lnB_st[w], E=E_st[w]),
                        lower=c(lnB0=-Inf, E=0),
                        upper=c(lnB0=Inf,  E=30),
                        weights=wgts,
                        control=list(minFactor=1 / 2^16, maxiter=1024)),
                  silent=TRUE)

    out = cbind(summary(arr_fit)$coefficients[, c('Estimate', 'Std. Error')],
                start=c(lnB_st[w], E_st[w]), vcov(arr_fit), AIC_c=aics_out[w],
                AIC=AIC(arr_fit), BIC=BIC(arr_fit), log_Lik=logLik(arr_fit))
  } else
  {
    out = matrix(as.numeric(NA), ncol=2 + 1 + 2 + 1, nrow=2)
  }

  par_names = c('lnB0', 'E')
  out = data.frame(par_names, out, stringsAsFactors=FALSE, row.names=NULL)
  colnames(out) = c('pars', 'est', 'se', 'start', 'lnB0', 'E', 'AIC_c', 'AIC',
                    'BIC', 'log_Lik')
  return(out)
}



School = function(lnB0, E, E_D, Th_K, Tref_K, tmp_K, k=BOLTZ)
{
  # Modified Schoolfield-Sharpe model (Schoolfield & Sharpe, 1981), used by Dan
  # and Gab in their ELett paper. Fitting is on ln trait scales.
  # ~~~ Parameters ~~~
  # lnB0   :numeric: Approximate ln rate at T_ref.
  # E      :numeric: Activation energy (eV).
  # E_D    :numeric: De-activation energy (eV).
  # Th_K   :numeric: T at which enzyme is 50% active and 50% high-temperature
  #         suppressed (in K). 
  # Tref_K :numeric: Reference temperature at which the rate is approximately
  #         ln_const (in K).
  # tmp_K  :numeric: temperatures (in K), as value or vector.
  # k      :numeric: Boltzmann constant.

  return(lnB0 - E/k * (1/tmp_K - 1/Tref_K) - 
         log(1 + exp((E_D/k) * (1/Th_K - 1/tmp_K))))
}



SchoolFit = function(dat, trait_name, temp_name, unc_name=NULL, Tref_K,
                     n_rand=1, k=BOLTZ)
{
  # Fit the Schoolfield model, and return NAs if no convergence. Optionally, use
  # randomised starting values, to improve chances of finding a fit. 
  # ~~~ Parameters ~~~
  # dat        : data.frame or data.table - including temperature and trait data.
  # trait_name : string - name of column that holds trait values.
  # temp_name  : string - name of column that holds temperature data (in K).
  # unc_name   : string - name of column with the measure of uncertainty in a
  #              point. The inverse of this will be used as weighting in the
  #              fitting. Otherwise, leave as NULL, and equal weighting will be
  #              used.
  # Tref_K     : reference temperature (in K).
  # n_rand     : numeric - number of random starting values used (if = 1, then
  #              no random starting values are used).
  # k          : Boltzmann constant.
  # ~~~ Output ~~~
  # Matrix with:
  # - Point estimates.
  # - Standard errors.
  # - Starting values used in the best fit.
  # - Variance covarience (4x4) matrix.
  # - AICc, AIC, BIC, and log likelihood values for teh fit.

  library('minpack.lm')
  library('sme')
  library('truncnorm')

  trait  = dat[[trait_name]]
  temp_K = dat[[temp_name]]

  if (length(unc_name) > 0)
  {
    wgts = 1 / dat[[unc_name]]
  } else
  {
    l    = length(trait)
    wgts = rep(1 / l, l)
  }

  if (min(temp_K, na.rm=TRUE) < 200)
  {
    stop('Error in SchoolFit : make sure temperatures are in K')
  }

  # Estimate Th as being approximately T_peak.
  Th_st  = .GetTpk(tmp_K=temp_K, rate=trait)
  # This is two values: the slope, and the standard error around slope.
  E_st   = as.numeric(.GetE(tmp_K=temp_K, rate=trait, Tp_K=Th_st))
  lnB_st = .GetLnB0(tmp_K=temp_K, rate=trait, Tref_K=Tref_K)
  B_st   = exp(lnB_st)

  # Create randomised starting points.
  E_st_pe = E_st[1]  # Slope value point estimate.
  Th_st   = c(Th_st, rnorm(n_rand-1, mean=Th_st, sd=15))

  if (n_rand > 1)
  {
    # We need truncated normal to ensure we don't get negative values of E.
    E_st = c(E_st_pe, rtruncnorm(n_rand-1, a=0, b=Inf, mean=E_st[1], sd=2 * E_st[2]))
    # Randomise on linear scale. Again, we don't want negative rates.
    lnB_st = c(lnB_st, log(rtruncnorm(n_rand-1, a=0, b=Inf, mean=B_st, sd=B_st / 2)))
  } else
  {
    E_st = E_st_pe
  }

  # Select the best model using AICc. Many of these turn out to be similar.
  aics_out = rep(NA, n_rand)

  for (i in 1:n_rand)
  {
    school_fit = try(nlsLM(log(trait) ~ School(lnB0, E, E_D, Th_K, 
                                               Tref_K=Tref_K, tmp_K=temp_K),
                           start=c(lnB0=lnB_st[i], E=E_st[i], E_D=4*E_st[i],
                                   Th_K=Th_st[i]),
                           lower=c(lnB0=-Inf, E=0,  E_D=0,  Th_K=250),
                           upper=c(lnB0=Inf,  E=30, E_D=50, Th_K=350),
                           weights=wgts,
                           control=list(minFactor=1 / 2^16, maxiter=1024)),
                     silent=TRUE)

    if (class(school_fit) != 'try-error')
    {
      aics_out[i] = AICc(school_fit)
    }
  }

  w = which.min(aics_out)

  if (length(w) > 0 & n_rand == 1)
  {
    out = cbind(summary(school_fit)$coefficients[, c('Estimate', 'Std. Error')],
                start=c(lnB_st, E_st, 4*E_st, Th_st), vcov(school_fit),
                AIC_c=aics_out, AIC=AIC(school_fit), BIC=BIC(school_fit),
                log_Lik=logLik(school_fit))
  } else if (length(w) > 0 & n_rand > 1)
  {
    school_fit = try(nlsLM(log(trait) ~ School(lnB0, E, E_D, Th_K, 
                                               Tref_K=Tref_K, tmp_K=temp_K),
                           start=c(lnB0=lnB_st[w], E=E_st[w], E_D=4*E_st[w],
                                   Th_K=Th_st[w]),
                           lower=c(lnB0=-Inf, E=0,  E_D=0,  Th_K=250),
                           upper=c(lnB0=Inf,  E=30, E_D=50, Th_K=350),
                           weights=wgts,
                           control=list(minFactor=1 / 2^16, maxiter=1024)),
                     silent=TRUE)
    out = cbind(summary(school_fit)$coefficients[, c('Estimate', 'Std. Error')],
                start=c(lnB_st[w], E_st[w], 4*E_st[w], Th_st[w]),
                vcov(school_fit), AIC_c=aics_out[w], AIC=AIC(school_fit),
                BIC=BIC(school_fit), log_Lik=logLik(school_fit))
  } else 
  {
    out = matrix(as.numeric(NA), ncol=2 + 1 + 4 + 1, nrow=4)
  }

  par_names = c('lnB0', 'E', 'E_D', 'Th_K')
  out = data.frame(par_names, out, stringsAsFactors=FALSE, row.names=NULL)
  colnames(out) = c('pars', 'est', 'se', 'start', 'lnB0', 'E', 'E_D', 'Th_K',
                    'AIC_c', 'AIC', 'BIC', 'log_Lik')
  return(out)
}



Tpeak = function(E, E_D, Th_K, k=BOLTZ)
{
  # Estimate the temperature at which a trait peaks from the Schoolfield
  # parameters.
  # ~~~ Parameters ~~~
  # E    :numeric: Activation energy (eV).
  # E_D  :numeric: De-activation energy (eV).
  # Th_K :numeric: T at which enzyme is 50% active and 50% high-temperature
  #       suppressed (in K). 
  # k    :numeric: Boltzmann constant.

  return(E_D * Th_K / (E_D + k * Th_K * log((E_D - E) / E)))
}



School_Tpk = function(lnB0, E, E_D, Tp_K, Tref_K, tmp_K, k=BOLTZ)
{
  # Modified Schoolfield-Sharpe model, in terms of Tpk instead of Th. 
  # ~~~ Parameters ~~~
  # lnB0     :numeric: Approximate ln rate at T_ref.
  # E        :numeric: Activation energy (eV).
  # E_D      :numeric: De-activation energy (eV).
  # Tp_K     :numeric: T at which rate is maximised (in K).
  # Tref_K   :numeric: Reference temperature at which the rate is approximately
  #           ln_const (in K).
  # tmp_K    :numeric: temperatures (in K), as value or vector.
  # k        :numeric: Boltzmann constant.

  return(lnB0 - E/k * (1/tmp_K - 1/Tref_K) - 
         log(1 + (E / (E_D-E)) * exp((E_D/k * (1/Tp_K - 1/tmp_K)))))
}


SchoolTpkFit = function(dat, trait_name, temp_name, unc_name=NULL, Tref_K,
                        n_rand=1, k=BOLTZ)
{
  # Fit the Tpk version of the Schoolfield model, and return NAs if no
  # convergence. Optionally, use randomised starting values, to improve chances
  # of finding a fit. 
  # ~~~ Parameters ~~~
  # dat        :data.frame: including temperature and trait data.
  # trait_name :string: name of column that holds trait values.
  # temp_name  :string: name of column that holds temperature data (in K).
  # unc_name   :string: name of column with the measure of uncertainty in a
  #             point. The inverse of this will be used as weighting in the
  #             fitting. Otherwise, leave as NULL, and equal weighting will be
  #             used.
  # Tref_K     :numeric: reference temperature (in K).
  # n_rand     :numeric: number of random starting values used (if = 1, then
  #             no random starting values are used).
  # k          :numeric: Boltzmann constant.
  # ~~~ Output ~~~
  # Matrix with:
  # - Point estimates.
  # - Standard errors.
  # - Starting values used in the best fit.
  # - Variance covarience (4x4) matrix.
  # - AICc, AIC, BIC, and log likelihood values for teh fit.

  library('minpack.lm')
  library('sme')
  library('truncnorm')

  trait  = dat[[trait_name]]
  temp_K = dat[[temp_name]]

  if (length(unc_name) > 0)
  {
    wgts = 1 / dat[[unc_name]]
  } else
  {
    l    = length(trait)
    wgts = rep(1 / l, l)
  }

  if (min(temp_K, na.rm=TRUE) < 200)
  {
    stop('Error in SchoolFit : make sure temperatures are in K')
  }

  # Estimate Th as being approximately T.peak.
  Tp_st  = .GetTpk(tmp_K=temp_K, rate=trait)
  # This is two values: the slope, and the standard error around slope.
  E_st   = as.numeric(.GetE(tmp_K=temp_K, rate=trait, Tp_K=Tp_st))
  lnB_st = .GetLnB0(tmp_K=temp_K, rate=trait, Tref_K=Tref_K)
  B_st   = exp(lnB_st)

  # Create randomised starting points.
  E_st_pe = E_st[1]  # Slope value point estimate.
  Tp_st   = c(Tp_st, rnorm(n_rand-1, mean=Tp_st, sd=15))

  if (n_rand > 1)
  {
    # We need truncated normal to ensure we don't get negative values of E.
    E_st = c(E_st_pe, rtruncnorm(n_rand-1, a=0, b=Inf, mean=E_st[1], sd=2 * E_st[2]))
    # Randomise on linear scale. Again, we don't want negative rates.
    lnB_st = c(lnB_st, log(rtruncnorm(n_rand-1, a=0, b=Inf, mean=B_st, sd=B_st / 2)))
  } else
  {
    E_st = E_st_pe
  }

  # Select the best model using AICc. Many of these turn out to be similar.
  aics_out = rep(NA, n_rand)

  for (i in 1:n_rand)
  {
    school_fit = try(nlsLM(log(trait) ~ School_Tpk(lnB0, E, E_D, Tp_K, 
                                                   Tref_K=Tref_K, tmp_K=temp_K),
                           start=c(lnB0=lnB_st[i], E=E_st[i], E_D=4*E_st[i],
                                   Tp_K=Tp_st[i]),
                           lower=c(lnB0=-Inf, E=0,  E_D=0,  Tp_K=250),
                           upper=c(lnB0=Inf,  E=30, E_D=50, Tp_K=350),
                           weights=wgts,
                           control=list(minFactor=1 / 2^16, maxiter=1024)),
                     silent=TRUE)

    if (class(school_fit) != 'try-error')
    {
      aics_out[i] = AICc(school_fit)
    }
  }

  w = which.min(aics_out)

  if (length(w) > 0 & n_rand == 1)
  {
    out = cbind(summary(school_fit)$coefficients[, c('Estimate', 'Std. Error')],
                start=c(lnB_st, E_st, 4*E_st, Tp_st), vcov(school_fit), 
                AIC_c=aics_out, BIC=BIC(school_fit), log_Lik=logLik(school_fit))
  } else if (length(w) > 0 & n_rand > 1)
  {
    school_fit = try(nlsLM(log(trait) ~ School(lnB0, E, E_D, Tp_K, 
                                               Tref_K=Tref_K, tmp_K=temp_K),
                           start=c(lnB0=lnB_st[w], E=E_st[w], E_D=4*E_st[w],
                                   Tp_K=Tp_st[w]),
                           lower=c(lnB0=-Inf, E=0,  E_D=0,  Tp_K=250),
                           upper=c(lnB0=Inf,  E=30, E_D=50, Tp_K=350),
                           weights=wgts,
                           control=list(minFactor=1 / 2^16, maxiter=1024)),
                     silent=TRUE)
    out = cbind(summary(school_fit)$coefficients[, c('Estimate', 'Std. Error')],
                start=c(lnB_st[w], E_st[w], 4*E_st[w], Tp_st[w]),
                vcov(school_fit), AIC_c=aics_out[w], BIC=BIC(school_fit),
                log_Lik=logLik(school_fit))
  } else
  {
    out = matrix(as.numeric(NA), ncol=2 + 1 + 4 + 1, nrow=4)
  }

  par_names = c('lnB0', 'E', 'E_D', 'Tp_K')
  out = data.frame(par_names, out, stringsAsFactors=FALSE, row.names=NULL)
  colnames(out) = c('pars', 'est', 'se', 'start', 'lnB0', 'E', 'E_D', 'Tp_K',
                    'AIC_c', 'BIC', 'log_Lik')
  return(out)
}



SchoolLowIn = function(lnB0, E, E_l_D, E_h_D, Th_l_K, Th_h_K, Tref_K, tmp_K,
                       k=BOLTZ)
{
  # Modified Schoolfield-Sharpe model (Schoolfield & Sharpe, 1981), including
  # low temperature inactivation.
  # ~~~ Parameters ~~~
  # lnB0   :numeric: Approximate ln rate at T_ref.
  # E      :numeric: Activation energy (eV).
  # E_h_D  :numeric: High temperature de-activation energy (eV).
  # E_l_D  :numeric: Low temperature de-activation energy (eV).
  # Th_h_K :numeric: Low temperature Th_K (K).
  # Th_l_K :numeric: Low temperature Th_K (K).
  # Tref_K :numeric: Reference temperature at which the rate is ln_const (in K).
  # tmp_K  :numeric: temperatures (in K), as value or vector.
  # k      :numeric: Boltzmann constant.

  return(lnB0 - E/k * (1/tmp_K - 1/Tref_K) - 
         log(1 + exp((E_h_D/k) * (1/Th_h_K - 1/tmp_K)) + 
             exp((E_l_D/k) * (1/Th_l_K - 1/tmp_K))))
}



SchoolLowInFit = function(dat, trait_name, temp_name, Tref_K, wgts_name=NULL,
                          n_rand=1, k=BOLTZ)
{
  # NOTE: I have not used this in ages - probably does not work.
  #
  # Fit the Schoolfield model, and return NAs if no convergence.
  # ~~~ Parameters ~~~
  # dat        :data.frame: including temperature and trait data.
  # trait_name :string: name of column that holds trait values.
  # temp_name  :string: name of column that holds temperature data (in K).
  # wgts_name  :string: name of column that holds the weights to be used. The
  #             inverse of the weights will be used in the fitting. Otherwise,
  #             leave as NULL, and equal weighting will be used.
  # Tref_K     :numeric: reference temperature (in K).
  # n_rand     :numeric: number of random starting values used (if = 1, then
  #             no random starting values are used).
  # k          :numeric: Boltzmann constant.

  library('minpack.lm')
  library('sme')
  library('truncnorm')

  trait  = dat[[trait_name]]
  temp_K = dat[[temp_name]]

  if (min(temp_K, na.rm=TRUE) < 200)
  {
    stop('Error in SchoolLowInFit : make sure temperatures are in K')
  }

  if (length(wgts_name) > 0)
  {
    wgts = 1 / dat[[wgts_name]]
  } else
  {
    l    = length(trait)
    wgts = rep(1 / l, l)
  }

  # Estimate Th as being approximately T_peak.
  Th_st  = .GetTpk(tmp_K=temp_K, rate=trait)
  E_st   = as.numeric(.GetE(tmp_K=temp_K, rate=trait, Tp_K=Th_st))
  lnB_st = .GetLnB0(tmp_K=temp_K, rate=trait, Tref_K=Tref_K)
  B_st   = exp(lnB_st)

  E_st_pe = E_st[1]  # Slope value point estimate.
  Th_st   = c(Th_st, rnorm(n_rand-1, mean=Th_st, sd=15))

  if (n_rand > 1)
  {
    E_st = c(E_st_pe, rtruncnorm(n_rand-1, a=0, b=Inf, mean=E_st[1], sd=2 * E_st[2]))
    lnB_st = c(lnB_st, log(rtruncnorm(n_rand-1, a=0, b=Inf, mean=B_st, sd=B_st / 2)))
  } else
  {
    E_st = E_st_pe
  }

  aics_out = rep(NA, n_rand)

  for (i in 1:n_rand)
  {
    s_fit = try(nlsLM(log(trait) ~ SchoolLowIn(lnB0, E, E_l_D, E_h_D, Th_l_K,
                                               Th_h_K, Tref_K=Tref_K,
                                               tmp_K=temp_K),
                      start=c(lnB0=lnB_st[i], E=E_st[i], E_h_D=4*E_st[i],
                              E_l_D=-3*E_st[i], Th_l_K=Th_st[i]-10,
                              Th_h_K=Th_st[i]),
                      lower=c(lnB0=-Inf, E=0, E_h_D=0, E_l_D=-Inf, Th_l_K=0,
                              Th_h_K=0),
                      upper=c(lnB0=Inf, E=Inf, E_h_D=Inf, E_l_D=0, Th_l_K=Inf,
                              Th_h_K=Inf),
                      weights=wgts,
                      control=list(minFactor=1 / 2^16, maxiter=1024)),
                silent=TRUE)

    if (class(s_fit) != 'try-error')
    {
      aics_out[i] = AICc(s_fit)
    }
  }

  w = which.min(aics_out)

  if (length(w) > 0 & n_rand == 1)
  {
    out = as.list(summary(s_fit)$coefficients[, c('Estimate', 'Std. Error')])
    out = c(lnB_st[w], E_st[w], 4*E_st[w], -3*E_st[w], Th_st[w]-10, Th_st[w],
            out, AIC(s_fit), AICc(s_fit), length(trait))
    names(out) = c(paste0(rep(names(coefficients(s_fit)), 2), 
                          rep(c('_st', '', '_se'), each=6)), 'aic', 'aic_c', 'N')
    return(out)
  } else if (length(w) > 0 & n_rand > 1)
  {
    s_fit = try(nlsLM(log(trait) ~ SchoolLowIn(lnB0, E, E_l_D, E_h_D, Th_l_K,
                                               Th_h_K, Tref_K=Tref_K,
                                               tmp_K=temp_K),
                      start=c(lnB0=lnB_st[w], E=E_st[w], E_h_D=4*E_st[w],
                              E_l_D=-3*E_st[w], Th_l_K=Th_st[w]-10,
                              Th_h_K=Th_st[w]),
                      lower=c(lnB0=-Inf, E=0, E_h_D=0, E_l_D=-Inf, Th_l_K=0,
                              Th_h_K=0),
                      upper=c(lnB0=Inf, E=Inf, E_h_D=Inf, E_l_D=0, Th_l_K=Inf,
                              Th_h_K=Inf),
                      weights=wgts,
                      control=list(minFactor=1 / 2^16, maxiter=1024)),
                silent=TRUE)

    out = as.list(summary(s_fit)$coefficients[, c('Estimate', 'Std. Error')])
    out = c(lnB_st[w], E_st[w], 4*E_st[w], -3*E_st[w], Th_st[w]-10, Th_st[w],
            out, AIC(s_fit), AICc(s_fit), length(trait))
    names(out) = c(paste0(rep(names(coefficients(s_fit)), 2), 
                          rep(c('_st', '', '_se'), each=6)), 'aic', 'aic_c', 'N')
    return(out)
  } else
  {
    return(as.list(as.numeric(c(lnB0_st=NA, E_st=NA, E_h_D_st=NA, E_l_D_st=NA,
                                Th_l_K_st=NA, Th_h_K_st=NA,
                                lnB0=NA, E=NA, E_h_D=NA, E_l_D=NA, Th_l_K=NA,
                                Th_h_K=NA,
                                lnB0_se=NA, E_se=NA, E_h_D_se=NA, E_l_D_se=NA,
                                Th_l_K_se=NA, Th_h_K_se=NA,
                                aic=NA, aic_c=NA, N=NA))))
  }
}



