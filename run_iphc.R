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

## The function run.logbook is specifically setup to run the CPUE
## standardization and hook spacing effect estimation on a data set. Here
## are the key inputs that control the model.

## form: 1=smoother effect; 2= Hamley and Skud formula; 3= constant (no effect)
## mdoel: 'NS'=not spatial; 'S'=spatial' ; 'ST'=full spatiotemporal
## vessel effect: TRUE=estimate vessel random effects; FALSE=ignore effect
## n_knots: controls the resolution of the spatial grid (~ # of cells)

## Loop through each regarea and get CPUE from full model to compare with
## survey
regareas <- '3A'
regareas <- c('2A', '2B', '2C', '3A', '3B', '4A', '4B')
knots <- 50
fits.areas <- lapply(regareas, function(x){
  d <- droplevels(subset(data.full, regcde==x))
  xx <- run.logbook(dat=d, n_knots=knots, model='NS', form=2,
                    vessel_effect=FALSE)
  return(xx)
})
## saveRDS(fits.areas, file='results/fits.areas.RDS')
## fits.areas <- readRDS('results/fits.areas.RDS')

## Load the results in and prepare data for plotting
survey <- readRDS('data/survey.RDS')
cpue.areas <- ldply(seq_along(regareas), function(x)
  cbind(data='logbook', regarea=regareas[x], fits.areas[[x]]$sd.density))
cpue.areas$se <- cpue.areas$sd
cpue.areas <- ddply(cpue.areas, .(regarea), mutate,
      cpue=value/mean(value),
      lwr=lwr/mean(value),
      upr=upr/mean(value))
cpue.areas <- cpue.areas[, names(survey)]
cpue.df <- rbind.fill(survey, cpue.areas)


g <- ggplot(cpue.df, aes(year, cpue, ymin=lwr, ymax=upr, group=data, fill=data)) +
  geom_ribbon(alpha=.5)+facet_wrap('regarea', scales='free_y')
g
## ggsave('plots/cpue_trends_vs_survey.png', g, width=ggwidth, height=ggheight)
