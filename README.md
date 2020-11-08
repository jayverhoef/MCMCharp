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

MCMC sampling for simple harp seal model with a 2-age Leslie matrix model where parameters are stochastic from year to year, but each parameter has a common prior across years.  There are also priors on starting values for adults and pups, and a prior that controls density-dependent decay in fecundity.

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

Run R Scripts
-------------

*Running MCMC code*

The main file that runs the MCMC code is the function MCMCharp().  To see a list of options type

```
help(MCMCharp)
```
To accept all defaults, and set a random number seed for repeatability, try

```
     W = MCMCharp(
       harvadu_data = catch_data$adu, 
       harvpup_data = catch_data$pup,
       pupcount_data = pup_production, 
       set_seed = 1001
     )
```

The results of this run are already stored, so if you don't want to wait, you can simply type

```
data(W)
```

*Some Graphics*



-------------
##### Disclaimer

<sub>This repository is a scientific product and is not official communication of the Alaska Fisheries Science Center, the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. AFSC-MML has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.</sub>
