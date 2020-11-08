library(MCMCharp)
data(catch_data)
data(pup_production)

#use all defaults
W = MCMCharp(
  harvadu_data = catch_data$adu, 
  harvpup_data = catch_data$pup,
  pupcount_data = pup_production, 
  set_seed = 1001
)

path1 = '/media/jay/data/desktop_data/ActiveService/NAMMCO/MCMCharp_package/'
save(W, file = paste0(path1,'MCMCharp/data/W.rda'))

