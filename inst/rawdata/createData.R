path1 = '/media/jay/data/desktop_data/ActiveService/NAMMCO/MCMCharp_package/'
catch_data = read.csv(paste0(path1,'MCMCharp/inst/rawdata/catch_data.csv'), 
  header = FALSE, col.names = c('year','pup','adu'))
pup_production = read.csv(paste0(path1,'MCMCharp/inst/rawdata/pup_production.csv'),
  header = FALSE, col.names = c('year','est','cv'))
fecundity_data = read.csv(paste0(path1,'MCMCharp/inst/rawdata/fecundity_data.csv'),
  header = FALSE, col.names = c('year','est','sd'))
pup_production$se = pup_production$cv*pup_production$est
catch_data$i = catch_data$year - 1944
fecundity_data$i = fecundity_data$year - 1944
pup_production$i = pup_production$year - 1944
save(catch_data, file = paste0(path1,'MCMCharp/data/catch_data.rda'))
save(fecundity_data, file = paste0(path1,'MCMCharp/data/fecundity_data.rda'))
save(pup_production, file = paste0(path1,'MCMCharp/data/pup_production.rda'))


