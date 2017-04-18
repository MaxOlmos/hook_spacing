### ------------------------------------------------------------
## Some quick code to look at other areas that the IPHC was interested
## in. NOT IN MANUSCRIPT.
source('startup.R')
## Load data for all regulatory areas, to be subsetted later
data.full <- readRDS(file='data/data.RDS')

Version <- "models/spatiotemporal_cpue_spacing"
clean.TMB.files(Version)
compile( paste0(Version,".cpp"))
dyn.load( dynlib(Version))

## Loop through each regarea and get CPUE from full model to compare with
## survey
regareas <- c('2A', '2B', '2C', '3A', '3B', '4A', '4B')
knots <- 500
fits.areas <- lapply(regareas, function(x){
  d <- droplevels(subset(data.full, regcde==x))
  xx <- run.logbook(dat=d, n_knots=knots, model='ST', form=2,
                    vessel_effect=TRUE)
  return(xx)
})
saveRDS(fits.areas, file='results/fits.areas.RDS')

