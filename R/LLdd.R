#-------------------------------------------------------------------------------
#
#          LL
#
#-------------------------------------------------------------------------------

#' loglikelihood for model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param N_pup1 population value for pups in first year
#' @param N_adu1 population value for adults in first year
#' @param delta Leslie matrix parameter for adult survival
#' @param phi Leslie matrix parameter for adult fecundity (both sexes, until survey time)
#' @param kappa1  Leslie matrix parameter for pup survival after surveys
#' @param rho  parameter for density dependent decay in fecundity
#' @param pdelta_mu mean of normal distribution for prior for logit(delta) 
#' @param pdelta_sd standard deviation of normal distribution for prior for logit(delta) 
#' @param pphi_mu mean of normal distribution for prior for logit(phi) 
#' @param pphi_sd standard deviation of normal distribution for prior for logit(phi) 
#' @param pkappa_mu mean of normal distribution for prior for logit(kappa)
#' @param pkappa_sd standard deviation of normal distribution for prior for logit(kappa) 
#' @param prho_mu mean of normal distribution for prior for log(rho)
#' @param prho_sd standard deviation of normal distribution for prior for log(rho)
#' @param pN_pup1_mu mean of normal distribution for prior for N_pup1
#' @param pN_pup1_sd standard devitation of normal distribution for prior for N_pup1
#' @param pN_adu1_mu mean of normal distribution for prior for N_adu1
#' @param pN_adu1_sd standard devitation of normal distribution for prior for N_adu1
#' @param pp  data set with pup production
#' @param H_adu vector of adult harvest for each year
#' @param H_pup vector of pup harvest for each year
#' @param nyrs  number of years (one more than length of harvest vectors)
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export
LLdd <- function(N_pup1, N_adu1, delta, phi, kappa1, rho,
        pdelta_mu, pdelta_sd, pphi_mu, pphi_sd, pkappa_mu, pkappa_sd, 
        prho_mu, prho_sd, pN_pup1_mu, pN_pup1_sd, pN_adu1_mu, pN_adu1_sd, 
        pp, H_adu, H_pup, nyrs)
{
  N_adu <- rep(NA, times = nyrs)
  N_pup <- N_adu
  N_adu[1] = N_adu1
  N_pup[1] = N_pup1
  delta = expit(delta)
  phi = expit(phi)
  kappa1 = expit(kappa1)
  rho = exp(rho)
  for(i in 2:nyrs) {
    N_adu[i] <- delta[i-1]*(N_adu[i-1]-H_adu[i-1]) + 
      kappa1[i-1]*(N_pup[i-1] - H_pup[i-1])
    N_pup[i] <- exp(-rho*(N_adu[i-1]+N_pup[i-1])/100000)*phi[i-1]*N_adu[i-1]
  }
  if(N_adu1 <= 0 | any((N_adu[2:nyrs]-H_adu) <= 0) | 
    N_pup1 <= 0 | any((N_pup[2:nyrs] - H_pup) <= 0)) return(-1e+32)
  sum(dnorm(pp[,"est"], N_pup[pp$i], pp[,"se"], log = TRUE)) +
      dnorm(N_pup1, pN_pup1_mu, pN_pup1_sd, log = TRUE) +  
      dnorm(N_adu1, pN_adu1_mu, pN_adu1_sd, log = TRUE) +
      sum(dnorm(logit(delta), pdelta_mu, pdelta_sd, log = TRUE)) +
      sum(dnorm(logit(phi), pphi_mu, pphi_sd, log = TRUE)) +
      sum(dnorm(logit(kappa1), pkappa_mu, pkappa_sd, log = TRUE)) +
      dnorm(log(rho), prho_mu, prho_sd, log = TRUE)
}
