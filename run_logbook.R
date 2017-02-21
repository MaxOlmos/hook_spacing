### ------------------------------------------------------------
## Step 1: read in and prep the data for this model
source('startup.R')
## Load data for all regulatory areas, to be subsetted later
data.full <- readRDS(file='data/data.RDS')

## data.temp <- ddply(data, .(geartype, year), summarize,
##                        total.catch=sum(catch))
## data.summarized <- ddply(data.temp, .(year), mutate, pct.catch=total.catch/sum(total.catch))
## saveRDS(data.summarized, file='results/data.summarized.RDS')

## data.simulated <- data
## data.simulated$catch <- with(data, hooks*exp(-.1+rnorm(n=nrow(data), mean=0, sd=.1)))
## data <- data.simulated
Version <- "models/spatiotemporal_cpue_spacing"
clean.TMB.files(Version)
compile( paste0(Version,".cpp"))
dyn.load( dynlib(Version))

### ------------------------------------------------------------
## Step 2. Model exploration. Models= no spatial effect (NS), spatiaggl model (S)
## and full spatio-temporal (ST). orm=1 implies a random walk on hook
## spacing, form=2 is the parametric HS model.

### Explore effects of key dimensions.
## Spacing vs model
knots <- 100
d <- droplevels(subset(data.full, regcde=='3A'))
vessel <- FALSE
fit1 <- run.logbook(d, n_knots=knots, model='NS', form=1, vessel=vessel)
fit2 <- run.logbook(d, n_knots=knots, model='NS', form=2, vessel=vessel)
fit3 <- run.logbook(d, n_knots=knots, model='NS', form=3, vessel=vessel)
fit4 <- run.logbook(d, n_knots=knots, model='ST', form=1, vessel=vessel)
fit5 <- run.logbook(d, n_knots=knots, model='ST', form=2, vessel=vessel)
fit6 <- run.logbook(d, n_knots=knots, model='ST', form=3, vessel=vessel)
fits.all <- list(fit1, fit2, fit3, fit4, fit5, fit6)
saveRDS(fits.all, file='results/fits.all.RDS')

## Loop through each regarea and get CPUE from full model to compare with
## survey
regareas <- c('2A', '2B', '2C', '3A', '3B', '4A', '4B')
knots <- 50
fits.areas <- lapply(regareas, function(x){
  d <- droplevels(subset(data.full, regcde==x))
  xx <- run.logbook(data=d, n_knots=knots, model='ST', form=2,
                    vessel=TRUE)
  return(xx)
})
saveRDS(fits.areas, file='results/fits.areas.RDS')

### Fit to some simulated data. Base if off the full model results from 3A.
fit <- readRDS('results/fits.areas.RDS')[[4]]
d <- droplevels(subset(data.full, regcde=='3A'))
fits.sim <- ldply(1:20, function(i)
  simulate.fit(i=i, d=d, fit=fit, knots=1000, model='ST'))
saveRDS(fits.sim, file='results/fits.sim.RDS')
ggplot(fits.sim, aes(year, rel.error, group=rep)) + geom_line() + facet_wrap('trend')



## Cleanup
dyn.unload( dynlib(Version))


### ------------------------------------------------------------
### OLD CODE

## ## Test code. Track parameter traces for plotting
## Inputs <- make.inputs(n_knots=50, model='S', form=1, likelihood=1)
## Obj <- MakeADATAun(data=Inputs$Data, parameters=Inputs$Params,
##                  random=Inputs$Random, map=Inputs$Map)
## Obj$env$beSilent()
## Opt <<- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
##               control=list(trace=0, eval.max=1e4, iter.max=1 ))

## xx <- ldply(1:200, function(i){
##               print(i)
## Opt <<- nlminb(start=Opt$par, objective=Obj$fn, gradient=Obj$gr,
##               control=list(trace=0, eval.max=1e4, iter.max=5 ))
## report.temp <- Obj$report();
## return(c(iteration=i, unlist(report.temp[c(2,6,7,8)]), Opt$par))
## })
## xx.long <- melt(xx, 'iteration')
## ggplot(xx.long, aes(iteration, value)) + geom_line() +
##   facet_wrap('variable', scales='free_y')

## ## Loop through each model, running and saving results.
## n_knots <- 50
## Opt.list <- Report.list <- SD.list <- list()
## for(m in c('NS', "S", "ST")){
##   print(paste0('Starting model: ', m))
##   Inputs <- make.inputs(n_knots=n_knots, model=m, likelihood=1)
##   Obj <- MakeADATAun( data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random,
##                    map=Inputs$Map)
##   trash <- Obj$env$beSilent()
##   start <- Sys.time()
##   temp <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
##                  control=list(trace=50, eval.max=1e4, iter.max=1e4))
##   Opt.list[[m]] <- nlminb( start=temp$par, objective=Obj$fn, gradient=Obj$gr,
##                      control=list(trace=50, eval.max=1e4, iter.max=1e4))
##   Opt.list[[m]][["final_diagnostics"]] <-
##     data.frame( "Name"=names(Obj$par), "final_gradient"=as.numeric(Obj$gr(Opt.list[[m]]$par)))
##   Report.list[[m]] <- test <- Obj$report()
##   Report.list[[m]]$time <- as.numeric(difftime(Sys.time(),start,
##                                                units='mins'))
##   SD.list[[m]] <- sdreport(Obj)
##   ##  make.model.plots(

## }
## saveRDS(Report.list, 'results/report_models.RDS')
## saveRDS(SD.list, 'results/sd_models.RDS')

###
## mapdata <- map_data("worldHires")
## mapdata <- mapdata[mapdata$region %in% c('USA', 'Canada'),]
## map.xlim <- c(-156, -137)
## map.ylim <- c(55,63)
## g <- ggplot()+geom_map(data=mapdata,map=mapdata, aes(x=long, y=lat, map_id=region)) +
##   theme(text=element_text(size=12))+
##     coord_cartesian(xlim=map.xlim, ylim=map.ylim)
## g2 <-
##   g+geom_point(data=data, aes(x=longitude, y=latitude, colour=statarea),
##                size=map.size, alpha=map.alpha)
## ggsave('plots/map_statareas.png', g2, width=map.width, height=map.height)
## png('plots/pairs.png', width=ggwidth, height=ggheight, units='in', res=500)
## pairs(data[, c('logcpue', 'hooksize', 'geartype', 'depth', 'month')])
## dev.off()


## plot(resids~spacing, data=temp2)

## ## Increase spatial resolution and see what happens
## knots <- c(100, 150, 200, 300, 500, 800, 1000)
## Opt.list <- Report.list <- SD.list <- list()
## m <- 'M3'
## for(k in knots){
##   print(paste0('Starting knots: ', k))
##   Inputs <- make.inputs(k, m)
##   Obj <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params,
##                    random=Inputs$Random, map=Inputs$Map)
##   trash <- Obj$env$beSilent()
##   start <- Sys.time()
##   temp <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
##                  control=list(trace=50, eval.max=1e4, iter.max=1e4))
##   Opt.list[[k]] <- nlminb( start=temp$par, objective=Obj$fn, gradient=Obj$gr,
##                      control=list(trace=50, eval.max=1e4, iter.max=1e4))
##   Opt.list[[k]][["final_diagnostics"]] <-
##     data.frame( "Name"=names(Obj$par),
##                "final_gradient"=as.numeric(Obj$gr(Opt.list[[k]]$par)))
##   Report.list[[k]] <- Obj$report()
##   Report.list[[k]]$time <-
##       as.numeric(difftime(Sys.time(),start, units='mins'))
##   SD.list[[k]] <- sdreport(Obj)
## }
## Results <- do.call(rbind, lapply(Report.list, function(x)
##   data.frame(x[c('intercept', 'beta_depth','SigmaE', 'Range', 'Sigma',
##                  'SigmaO', 'jnll', 'time')])))
## Results$knots <- knots[1:7]
## Results.long <- melt(Results, id.vars='knots')
## g <- ggplot(Results.long, aes(knots, (value), group=variable)) +
##   facet_wrap('variable', scales='free_y') + geom_line()
## g
## ggsave('plots/results_by_resolution.png', g, width=ggwidth, height=ggheight)


