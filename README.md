[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.1.1-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kotzeb0912)](https://cran.r-project.org/package=kotzeb0912) [![packageversion](https://img.shields.io/badge/Package%20version-1.0-orange.svg?style=flat-square)](commits/master)

[![Last-changedate](https://img.shields.io/badge/last%20change-2020--11--08-yellowgreen.svg)](/commits/master)

# MCMCharp 
## A Quick and Dirty Bayesian Leslie Matrix Model of Norway HarpEast Data 

#### Jay M. Ver Hoef<sup>a</sup>

#### <sup>a</sup>NOAA Fisheries (NMFS) Alaska Fisheries Science Center

As a scientific work, and in keeping with common scientific practicies, we kindly request that you cite our research project and applicable publications if you use our work(s) or data in your publications or presentations. Additionally, we strongly encourage and welcome collaboration to promote use of these data in the proper context and scope.  The publication is currently in revision:


Executive Summary
-----------------

MCMC sampling for a simple harp seal model with a 2-age Leslie matrix model where parameters are stochastic from year to year, but each parameter has a common prior across years.  There are also priors on starting values for adults and pups, and a prior that controls density-dependent decay in fecundity.

Installation
------------

Installation of this R data package is done through the `devtools::install_github()` function or by downloading the [source package from the latest release](https://github.com/jayverhoef/MCMCharp).

```
library("devtools")
install_github("jayverhoef/MCMCharp")
```

Examine the Example Data
------------------------

The data are embedded in the R package

```
library(MCMCharp)
data(catch_data)
data(pup_production)
print(catch_data)
print(pup_production)
```
in R.  The scripts below show how to import it and obtain results in the paper.

The data were created from CSV files stored here

system.file("rawdata/catch_data.csv", package = "MCMCharp")

system.file("rawdata/pup_production.csv", package = "MCMCharp")

using this script

system.file("rawdata/createData.R", package = "MCMCharp")

You can navigate to this script and run it in you favorite R environment.

A PDF Presentation on Model and Results
------------------------

A pdf presentation on the model can be found here

system.file("doc/NAMMCO_Harpeast.pdf", package = "MCMCharp")

The LATEX file that produced the presentation can be found here

system.file("doc/NAMMCO_Harpeast.tex", package = "MCMCharp")

and all figures can be found in this subdirectory

system.file("doc/figure", package = "MCMCharp")

with the R script that created the figures here

system.file("doc/figure/create_pdf_figures.R", package = "MCMCharp")

Run R Scripts
-------------

*Running MCMC code*

The main file that runs the MCMC code is the function MCMCharp().  To see a list of options type

```
help(MCMCharp)
```

The arguments have descriptions, and match the description of the code, given above.

To accept all defaults, and set a random number seed for repeatability, try

```
     W = MCMCharp(
       harvadu_data = catch_data$adu, 
       harvpup_data = catch_data$pup,
       pupcount_data = pup_production, 
       set_seed = 1001
     )
```

This code, contained in an R script, can be found here,

```
system.file("scripts/runMCMC.R", package = "MCMCharp")
```

The results of this run are already stored, so if you don't want to wait, you can simply type to go straight to graphics below.

```
data(W)
```


*Some Graphics*

Here are some trace plots for MCMC sampling.  More formal tests on convergence can be applied.  A single R script file with all of these plots can be found here

system.file("scripts/create_figures.R", package = "MCMCharp")


trace plot for N_adu[1]

```
plot(W$N_adu1, type = 'l')
```

trace plot for N_pup[1]

```
plot(W$N_pup1, type = 'l')
```

some trace plots for phi for several years (indices 1, 2, 59, 60)

```
plot(unlist(lapply(W$phi, function(x) x[[1]])), type = 'l')
plot(unlist(lapply(W$phi, function(x) x[[2]])), type = 'l')
plot(unlist(lapply(W$phi, function(x) x[[59]])), type = 'l')
plot(unlist(lapply(W$phi, function(x) x[[60]])), type = 'l')
```

some trace plots for delta for several years (indices 1, 2, 59, 60)

```
plot(unlist(lapply(W$delta, function(x) x[[1]])), type = 'l')
plot(unlist(lapply(W$delta, function(x) x[[2]])), type = 'l')
plot(unlist(lapply(W$delta, function(x) x[[59]])), type = 'l')
plot(unlist(lapply(W$delta, function(x) x[[60]])), type = 'l')
```

some trace plots for kappa for several years (indices 1, 2, 59, 60)

```
plot(unlist(lapply(W$kappa1, function(x) x[[1]])), type = 'l')
plot(unlist(lapply(W$kappa1, function(x) x[[2]])), type = 'l')
plot(unlist(lapply(W$kappa1, function(x) x[[59]])), type = 'l')
plot(unlist(lapply(W$kappa1, function(x) x[[60]])), type = 'l')
```

trace plot for rho

```
plot(W$rho, type = 'l')
```

fitted abundance and trajectory

```
  H_adu = catch_data$adu
  H_pup = catch_data$pup
  pp = pup_production
  plot((1:75) + 1944, 1:75, ylim = c(0, 2000000), type = 'n',
    xlab = 'Year', ylab = 'Population')
    
  nyrs = 75
  for(k in 1:1000) {
    N_adu <- rep(NA, times = nyrs)
    N_pup <- N_adu
    N_adu[1] <- W$N_adu1[k]
    N_pup[1] <- W$N_pup1[k]
    for(i in 2:nyrs) {
      N_adu[i] <- W$delta[[k]][i-1]*(N_adu[i-1]-H_adu[i-1]) + 
        W$kappa1[[k]][i-1]*(N_pup[i-1] - H_pup[i-1])
      N_pup[i] <- exp(-W$rho[k]*(N_adu[i-1]+N_pup[i-1])/100000)*
        W$phi[[k]][i-1]*N_adu[i-1]
    }
  lines((1:length(N_adu)) + 1944, N_adu, col = rgb(.5,.5,.5,.05))
  lines((1:length(N_adu)) + 1944, N_pup, col = rgb(.9,.1,.1,.03))
  }
  
  points(pp[,c('year','est')], pch = 19, col = 'brown', cex = 2)
  lines(catch_data$year, catch_data$adu, lwd = 3)
  lines(catch_data$year, catch_data$pup, col = 'orange', lwd = 3)
```

priors for the following figures

```
pdelta_mu = 2.2
pdelta_sd = .5
pphi_mu = -.4
pphi_sd = .5
pkappa_mu = 1.4
pkappa_sd = .5
prho_mu = -4.0
prho_sd = .5
pN_pup1_mu = 300000
pN_pup1_sd =  50000
pN_adu1_mu = 1000000
pN_adu1_sd =  100000
```

Priors and Posteriors for phi for several years

```
  par_orig = par(no.readonly = TRUE)
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[1]])))),
    main = '1946', xlab = 'logit(phi[1])', lwd = 2)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[2]])))),
    main = '1947', xlab = 'logit(phi[2])', lwd = 2)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[59]])))),
    main = '2003', xlab = 'logit(phi[59])', lwd = 2)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[60]])))),
    main = '2004', xlab = 'logit(phi[60])', lwd = 2)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  par(par_orig)
```

priors and posteriors for delta for several years

```
  par_orig = par(no.readonly = TRUE)
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[1]])))),
    main = '1946', xlab = 'logit(delta[1])', lwd = 2)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[2]])))),
    main = '1947', xlab = 'logit(delta[2])', lwd = 2)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[59]])))),
    main = '2003', xlab = 'logit(delta[59])', lwd = 2)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[60]])))),
    main = '2004', xlab = 'logit(delta[60])', lwd = 2)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  par(par_orig)
```

priors and posteriors for kappa for several years

```
  par_orig = par(no.readonly = TRUE)
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[1]])))),
    main = '1946', xlab = 'logit(kappa[1])', lwd = 2)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[2]])))),
    main = '1947', xlab = 'logit(kappa[2])', lwd = 2)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[59]])))),
    main = '2003', xlab = 'logit(kappa[59])', lwd = 2)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[60]])))),
    main = '2004', xlab = 'logit(kappa[60])', lwd = 2)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  par(par_orig)
```

prior and posterior for rho 

```
plot(density(log(W$rho)), xlim = c(-5.5, -1.5),
  xlab = 'log(rho)', lwd = 2)
lines((-30:30)/30 + prho_mu, 
  dnorm((-30:30)/30 + prho_mu, mean = prho_mu, sd = prho_sd), col = 'blue', lwd = 2)
```

priors and posteriors for first eigenvalue for several years, both with and without density dependence factor (without is intrinsic growth at very low population values)

```
  par_orig = par(no.readonly = TRUE)
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = unlist(lapply(W$phi, function(x) x[[1]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[1]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[1]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '1946, No Dens Dep',
    xlab = 'Posterior First Eigenvalue')
  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = exp(-W$rho[k]*(N_adu[i-1]+N_pup[i-1])/100000)*
      unlist(lapply(W$phi, function(x) x[[1]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[1]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[1]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '1946, Dens Dep',
    xlab = 'Posterior First Eigenvalue')   
  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = unlist(lapply(W$phi, function(x) x[[59]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[59]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[59]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '2003, No Dens Dep',
    xlab = 'Posterior First Eigenvalue')
  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = exp(-W$rho[k]*(N_adu[i-1]+N_pup[i-1])/100000)*
      unlist(lapply(W$phi, function(x) x[[59]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[59]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[59]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '2003, Dens Dep',
    xlab = 'Posterior First Eigenvalue')
  par(par_orig)
```

-------------
##### Disclaimer

<sub>This repository is a scientific product and is not official communication of the Alaska Fisheries Science Center, the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. AFSC-MML has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.</sub>
