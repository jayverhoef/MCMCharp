#-------------------------------------------------------------------------------
#
#          MCMC
#
#-------------------------------------------------------------------------------

#' MCMC sampling
#'
#' MCMC sampling for simple harp seal model with a 2-age Leslie matrix model where parameters are stochastic from year to year, but each parameter has a common prior across years.  There are also priors on starting values for adults and pups, and a prior that controls density-dependent decay in fecundity.
#'
#' @param burnin the number of MCMC samples for burnin. No MCMC samples will be saved until burnin has completed.
#' @param niter the number of MCMC iterations after burnin, when samples will be saved.
#' @param thin thinning of the MCMC chain.  This is the number of samples to skip over before saving a sample.
#' @param pdelta_mu mean of normal distribution for prior for logit(delta), where delta is the Leslie matrix parameter for adult survival. Default is 2.2
#' @param pdelta_sd standard deviation of normal distribution for prior for logit(delta). Default is 0.5
#' @param pphi_mu mean of normal distribution for prior for logit(phi), where phi is the Leslie matrix parameter for adult adult fecundity (both sexes, until survey time). Default is -0.4
#' @param pphi_sd standard deviation of normal distribution for prior for logit(phi). Default is 0.5
#' @param pkappa_mu mean of normal distribution for prior for logit(kappa), where kappa is the Leslie matrix parameter for pup survival after surveys. Default is 1.4
#' @param pkappa_sd standard deviation of normal distribution for prior for logit(kappa). Default is 0.5
#' @param prho_mu mean of normal distribution for prior for log(rho), where rho is a parameter for density dependent decay in fecundity. Default is -4.0
#' @param prho_sd standard deviation of normal distribution for prior for log(rho). Default is 0.5
#' @param pN_pup1_mu mean of normal distribution for prior for N_pup1, the population value for pups in first year. Default is 300000
#' @param pN_pup1_sd standard devitation of normal distribution for prior for N_pup1. Default is 50000
#' @param pN_adu1_mu mean of normal distribution for prior for N_adu1, the population value for adults in the first year. Default 1000000
#' @param pN_adu1_sd standard devitation of normal distribution for prior for N_adu1. Default is 100000
#' @param N_pup1_tune Tuning the Metropolis proposal, which is the standard deviation of a normal distribution centered on the current value for N_pup1. Default is 10000
#' @param N_adu1_tune Tuning the Metropolis proposal, which is the standard deviation of a normal distribution centered on the current value for N_adu1. Default is 10000
#' @param delta_tune Tuning the Metropolis proposal, which is the standard deviation of a normal distribution centered on the current value for logit(delta). Default is 0.05
#' @param phi_tune Tuning the Metropolis proposal, which is the standard deviation of a normal distribution centered on the current value for logit(phi). Default 0.05
#' @param kappa_tune Tuning the Metropolis proposal, which is the standard deviation of a normal distribution centered on the current value for logit(kappa). Default is 0.05
#' @param rho_tune Tuning the Metropolis proposal, which is the standard deviation of a normal distribution centered on the current value for log(rho). Default is 0.03
#' @param harvadu_data vector of adult harvest values for each year
#' @param harvpup_data vector of adult harvest values for each year
#' @param pupcount_data data set of pup counts
#' @param set_seed default If this is specified, the seed value for randomization will be set, so all results will be repeatable.  If NULL, the seed value is taken from the computer clock. Default is NULL.
#'
#' @return a list of the MCMC chain values (first six list items), as well as acceptance rates (last 6 list items).
#'
#' @author Jay Ver Hoef
#' @export
#' @examples
#' library(MCMCharp) 
#' data(catch_data)
#' data(pup_production)
#' #use all defaults
#' W = MCMCharp(
#'   harvadu_data = catch_data$adu, 
#'   harvpup_data = catch_data$pup,
#'   pupcount_data = pup_production, 
#'   set_seed = 1001
#' )
MCMCharp <- function(burnin = 10000, niter = 100000, thin = 100, 
  pdelta_mu = 2.2, pdelta_sd = 0.5, pphi_mu = -0.4, pphi_sd = 0.5,
  pkappa_mu = 1.4, pkappa_sd = 0.5, prho_mu = -4.0, prho_sd = 0.5,
  pN_pup1_mu = 300000, pN_pup1_sd =  50000, pN_adu1_mu = 1000000,
  pN_adu1_sd =  100000, N_pup1_tune = 10000, N_adu1_tune = 10000,
  delta_tune = 0.05, phi_tune = 0.05, kappa_tune = 0.05, rho_tune = 0.03,
  harvadu_data, harvpup_data, pupcount_data,
  set_seed = NULL)
{
  if(!is.null(set_seed)) set.seed(set_seed)
  # rename data sets to something shorter
  pp = pupcount_data
  H_adu = harvadu_data
  H_pup = harvpup_data
  nyrs = length(H_adu) + 1
  #create starting values near mean of priors
  delta = expit(rnorm(nyrs, pdelta_mu, pdelta_sd/100))
  phi = expit(rnorm(nyrs, pphi_mu, pphi_sd/100))
  kappa1 = expit(rnorm(nyrs, pkappa_mu, pkappa_sd/100))
  rho = exp(rnorm(1, prho_mu, prho_sd/100))
  N_adu1 = rnorm(1, pN_adu1_mu, pN_adu1_sd/100)
  N_pup1 = rnorm(1, pN_pup1_mu, pN_pup1_sd/100)

    for(k in 1:burnin) {
    # starting population for adults
    U <- log(runif(1))
    # make a proposal
    N_adu1.try <- rnorm(1, N_adu1, N_adu1_tune)
    # compute MH ratio on log scale (subtract)
    LLdif <- LLdd(N_pup1, N_adu1.try, logit(delta), logit(phi), logit(kappa1), 
        log(rho), pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    # keep if (log) ratio indicates
    if(LLdif > U) N_adu1 <- N_adu1.try

    # starting population for pups
    U <- log(runif(1))
    N_pup1.try <- rnorm(1, N_pup1, N_pup1_tune)
    LLdif <- LLdd(N_pup1.try, N_adu1, logit(delta), logit(phi), logit(kappa1), 
        log(rho), pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) N_pup1 <- N_pup1.try
    # delta
    U <- log(runif(1))
    delta.try <- rnorm(length(delta), logit(delta), delta_tune)
    LLdif <- LLdd(N_pup1, N_adu1, delta.try, logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) delta <- expit(delta.try)
    # phi
    U <- log(runif(1))
    phi.try <- rnorm(length(phi), logit(phi), phi_tune)
    LLdif <- LLdd(N_pup1, N_adu1, logit(delta), phi.try, logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) phi <- expit(phi.try)
    # kappa
    U <- log(runif(1))
    kappa1.try <- rnorm(length(kappa1), logit(kappa1), kappa_tune)
    LLdif <- LLdd(N_pup1, N_adu1, logit(delta), logit(phi), kappa1.try, log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) kappa1 <- expit(kappa1.try)
    # rho
    U <- log(runif(1))
    rho.try <- rnorm(1, log(rho), rho_tune)
    LLdif <- LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), rho.try,
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) rho <- exp(rho.try)

    if(k%%thin == 0) {
      cat("\r", "Burnin Number: ", k)
    }
  }
  
  W <- vector("list", 12)
  N_adu1_acc = 0
  N_pup1_acc = 0
  delta_acc = 0
  phi_acc = 0
  kappa_acc = 0
  rho_acc = 0

  for(k in 1:niter) {

    # starting population for adults
    U <- log(runif(1))
    N_adu1.try <- rnorm(1, N_adu1, N_adu1_tune)
    LLdif <- LLdd(N_pup1, N_adu1.try, logit(delta), logit(phi), logit(kappa1), 
        log(rho), pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) {
      N_adu1 <- N_adu1.try
      N_adu1_acc = N_adu1_acc + 1
    }
    # starting population for pups
    U <- log(runif(1))
    N_pup1.try <- rnorm(1, N_pup1, N_pup1_tune)
    LLdif <- LLdd(N_pup1.try, N_adu1, logit(delta), logit(phi), logit(kappa1), 
    log(rho), pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) {
      N_pup1 <- N_pup1.try
      N_pup1_acc = N_pup1_acc + 1
    }
    # delta
    U <- log(runif(1))
    delta.try <- rnorm(length(delta), logit(delta), delta_tune)
    LLdif <- LLdd(N_pup1, N_adu1, delta.try, logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) {
      delta <- expit(delta.try)
      delta_acc = delta_acc + 1
    }
    # phi
    U <- log(runif(1))
    phi.try <- rnorm(length(phi), logit(phi), phi_tune)
    LLdif <- LLdd(N_pup1, N_adu1, logit(delta), phi.try, logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) {
        phi <- expit(phi.try)
        phi_acc = phi_acc + 1
    }
    # kappa
    U <- log(runif(1))
    kappa1.try <- rnorm(length(kappa1), logit(kappa1), kappa_tune)
    LLdif <- LLdd(N_pup1, N_adu1, logit(delta), logit(phi), kappa1.try, log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) {
      kappa1 <- expit(kappa1.try)
      kappa_acc = kappa_acc + 1
    }
    # rho
    U <- log(runif(1))
    rho.try <- rnorm(1, log(rho), rho_tune)
    LLdif <- LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), rho.try,
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs) -
      LLdd(N_pup1, N_adu1, logit(delta), logit(phi), logit(kappa1), log(rho),
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
    if(LLdif > U) {
      rho <- exp(rho.try)
      rho_acc = rho_acc + 1
    }

    if(k%%thin == 0) {
      cat("\r", "Iteration Number: ", k)
      W[[1]] <- c(W[[1]],N_pup1)
      W[[2]] <- c(W[[2]],N_adu1)
      W[[3]] <- c(W[[3]],list(delta))
      W[[4]] <- c(W[[4]],list(phi))
      W[[5]] <- c(W[[5]],list(kappa1))
      W[[6]] <- c(W[[6]],rho)
    }
  }
  W[[7]] <- N_pup1_acc/niter
  W[[8]] <- N_adu1_acc/niter
  W[[9]] <- delta_acc/niter
  W[[10]] <- phi_acc/niter
  W[[11]] <- kappa_acc/niter
  W[[12]] <- rho_acc/niter

  cat("\n")
  names(W) <- c("N_pup1", "N_adu1", "delta", "phi", "kappa1", "rho",
  "N_pup1_acc", "N_adu1_acc", "delta_acc", "phi_acc", "kappa_acc", "rho_acc")
  W
}
