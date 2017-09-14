### ------------------------------------------------------------
## Step 1: read in and prep the data, and compile and link model
## If using real data, use this data set.
data.full <- readRDS(file='data/data.RDS')
d <- droplevels(subset(data.full, regcde=='3A'))
## ## Otherwise can simulate a set like this. See function for more control.
## data <- simulate.data(n_sites=10000, n_knots=200, beta=0.06)
## d <- data$dat                           # Simulated longline sets
## density <- data$density_t               # True relative abundance
## Compile and link TMB model
m <- "models/spatiotemporal_cpue_spacing"
## clean.TMB.files(m)
compile( paste0(m,".cpp"))
dyn.load( dynlib(m))

### ------------------------------------------------------------
## Step 2. Model exploration. Models= no spatial effect (NS), spatiaggl model (S)
## and full spatio-temporal (ST). form=1 implies a random walk on hook
## spacing, form=2 is the parametric HS model.

### Explore effcets of key dimensions.
## Spacing vs model
knots <- 2000
vessel <- TRUE
fit1 <- run.logbook(d, n_knots=knots, model='NS', form=1, vessel=vessel)
fit2 <- run.logbook(d, n_knots=knots, model='NS', form=2, vessel=vessel)
fit3 <- run.logbook(d, n_knots=knots, model='NS', form=3, vessel=vessel)
fit4 <- run.logbook(d, n_knots=knots, model='ST', form=1, vessel=vessel)
fit5 <- run.logbook(d, n_knots=knots, model='ST', form=2, vessel=vessel)
fit6 <- run.logbook(d, n_knots=knots, model='ST', form=3, vessel=vessel)
fits.all <- list(fit1, fit2, fit3, fit4, fit5, fit6)
table.runtime <- ldply(fits.all, summarize,
 model.name=model.name, form.name=form.name, runtime=round(runtime,2))
table.runtime
saveRDS(fits.all, file='results/fits.all.RDS')
## Cleanup
dyn.unload(dynlib(m))

## Test vessel effect
knots <- 10
v0 <- run.logbook(d, n_knots=knots, model='NS', form=2, vessel=FALSE)
v1 <- run.logbook(d, n_knots=knots, model='NS', form=2, vessel=TRUE)
v0$sd.par[30,2]
v1$sd.par[30,2]
## % reduction
(v0$sd.par[30,2]-v1$sd.par[30,2])/v1$sd.par[30,2]

### End of file
