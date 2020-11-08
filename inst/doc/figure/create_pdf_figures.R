library(MCMCharp)
data(catch_data)
data(pup_production)
data(W)

# path to local directory
# path = your_path
path = paste0('/media/jay/data/desktop_data/ActiveService/NAMMCO/',
  'MCMCharp_package/MCMCharp/inst/doc/figure/')

#fitted abundance and trajectory
pdf(file = paste0(path,'Pop_traj.pdf'))
  H_adu = catch_data$adu
  H_pup = catch_data$pup
  pp = pup_production
  par_orig = par(no.readonly = TRUE)
  par(mar = c(5,5,2,1))
  plot((1:75) + 1944, 1:75, ylim = c(0, 2000000), type = 'n',
    xlab = 'Year', ylab = 'Population', cex.lab = 2, cex.axis = 1.5)
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
  lines(catch_data$year, catch_data$adu, lwd = 4)
  lines(catch_data$year, catch_data$pup, col = 'orange', lwd = 4)
  par(par_orig)
dev.off()

#priors for the following figures
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

# Priors and Posteriors for phi for several years
pdf(file = paste0(path,'Post_phi.pdf'))
  par_orig = par(mar = c(5,5,5,1))
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[1]])))), cex.main = 2,
    main = '1946', xlab = 'logit(phi[1])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[2]])))), cex.main = 2,
    main = '1947', xlab = 'logit(phi[2])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[59]])))), cex.main = 2,
    main = '2003', xlab = 'logit(phi[59])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$phi, function(x) x[[60]])))), cex.main = 2,
    main = '2004', xlab = 'logit(phi[60])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pphi_mu, 
    dnorm((-30:30)/20 + pphi_mu, mean = pphi_mu, sd = pphi_sd), col = 'blue', lwd = 2)
  par(par_orig)
dev.off()

# Priors and Posteriors for delta for several years
pdf(file = paste0(path,'Post_delta.pdf'))
  par_orig = par(mar = c(5,5,5,1))
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[1]])))), cex.main = 2,
    main = '1946', xlab = 'logit(delta[1])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[2]])))), cex.main = 2,
    main = '1947', xlab = 'logit(delta[2])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[59]])))), cex.main = 2,
    main = '2003', xlab = 'logit(delta[59])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$delta, function(x) x[[60]])))), cex.main = 2,
    main = '2004', xlab = 'logit(delta[60])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pdelta_mu, 
    dnorm((-30:30)/20 + pdelta_mu, mean = pdelta_mu, sd = pdelta_sd), col = 'blue', lwd = 2)
  par(par_orig)
dev.off()

# Priors and Posteriors for kappa for several years
pdf(file = paste0(path,'Post_kappa.pdf'))
  par_orig = par(mar = c(5,5,5,1))
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[1]])))), cex.main = 2,
    main = '1946', xlab = 'logit(kappa[1])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[2]])))), cex.main = 2,
    main = '1947', xlab = 'logit(kappa[2])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[59]])))), cex.main = 2,
    main = '2003', xlab = 'logit(kappa[59])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  plot(density(logit(unlist(lapply(W$kappa1, function(x) x[[60]])))), cex.main = 2,
    main = '2004', xlab = 'logit(kappa[60])', lwd = 2, cex.lab = 2, cex.axis = 1.5)
  lines((-30:30)/20 + pkappa_mu, 
    dnorm((-30:30)/20 + pkappa_mu, mean = pkappa_mu, sd = pkappa_sd), col = 'blue', lwd = 2)
  par(par_orig)
dev.off()

# Prior and Posterior for rho 
pdf(file = paste0(path,'Post_rho.pdf'))
  par_orig = par(mar = c(5,5,5,1))
  plot(density(log(W$rho)), xlim = c(-5.5, -1.5),
    xlab = 'log(rho)', lwd = 3, cex.lab = 2, cex.axis = 2, cex.main = 1.5)
  lines((-30:30)/30 + prho_mu, 
    dnorm((-30:30)/30 + prho_mu, mean = prho_mu, sd = prho_sd), col = 'blue', 
    lwd = 3)
  par(par_orig)
dev.off()


# Priors and Posteriors for first eigenvalue for several years, both with and
# without density dependence factor (without is intrinsic growth at very low
# population values)

pdf(file = paste0(path,'Post_eigval.pdf'))
  par_orig = par(mar = c(5,5,5,1))
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = unlist(lapply(W$phi, function(x) x[[1]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[1]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[1]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '1946, No Dens Dep', cex.main = 1.5,
    xlab = 'Posterior First Eigenvalue', lwd = 3, cex.lab = 1.5, cex.axis = 1.3)

  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = exp(-W$rho[k]*(N_adu[i-1]+N_pup[i-1])/100000)*
      unlist(lapply(W$phi, function(x) x[[1]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[1]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[1]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '1946, Dens Dep', cex.main = 1.5,
    xlab = 'Posterior First Eigenvalue', lwd = 3, cex.lab = 1.5, cex.axis = 1.3)
   
  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = unlist(lapply(W$phi, function(x) x[[59]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[59]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[59]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '2003, No Dens Dep', cex.main = 1.5,
    xlab = 'Posterior First Eigenvalue', lwd = 3, cex.lab = 1.5, cex.axis = 1.3)

  evals = 1:1000
  for(k in 1:1000) {
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,2] = exp(-W$rho[k]*(N_adu[i-1]+N_pup[i-1])/100000)*
      unlist(lapply(W$phi, function(x) x[[59]]))[k]
    M[2,1] = unlist(lapply(W$kappa1, function(x) x[[59]]))[k]
    M[2,2] = unlist(lapply(W$delta, function(x) x[[59]]))[k]
    evals[k] = eigen(M)$values[1]
  }
  plot(density(evals), main = '2003, Dens Dep', cex.main = 1.5,
    xlab = 'Posterior First Eigenvalue', lwd = 3, cex.lab = 1.5, cex.axis = 1.3)
  par(par_orig)
dev.off()

