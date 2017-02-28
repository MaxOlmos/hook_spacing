### ------------------------------------------------------------
## Step 1: read in and prep the data, and compile and link model
source('startup.R')
## If using real data, use this data set.
data.full <- readRDS(file='data/data.RDS')
d <- droplevels(subset(data.full, regcde=='3A'))
## Otherwise can simulate a set like this. See function for more control.
data <- simulate.data(n_sites=10000, n_knots=200, beta=0.06)
d <- data$dat                           # Simulated longline sets
density <- data$density_t               # True relative abundance
## Compile and link TMB model
m <- "models/spatiotemporal_cpue_spacing"
## clean.TMB.files(m)
compile( paste0(m,".cpp"))
dyn.load( dynlib(m))

### ------------------------------------------------------------
## Step 2. Model exploration. Models= no spatial effect (NS), spatiaggl model (S)
## and full spatio-temporal (ST). orm=1 implies a random walk on hook
## spacing, form=2 is the parametric HS model.

### Explore effcets of key dimensions.
## Spacing vs model
knots <- 10
vessel <- FALSE
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

## Run test of resolution. Takes a very long time!
if(FALSE){
## Increase spatial resolution and see what happens
knots <- c(100, 150, 200, 300, 500, 800, 1000, 1200, 1500, 2000, 2500)
fits.res <- list()
for(i in 1:length(knots)){
  k <- knots[i]
  print(paste0('Starting knots: ', k))
  fits.res[[i]] <-
    run.logbook(d, n_knots=k, model='ST', form=2, vessel=FALSE)
  ## Might crash so save at each iteration
  saveRDS(fits.res, file='results/fits.res.RDS')
}
fits.res <- readRDS('results/fits.res.RDS')
res <- do.call(rbind, lapply(fits.res, function(x)
  cbind(knots=x$n_knots, x$sd.par[,1:3])))
res <- res[res$par %in% c('intercept', 'beta_depth','SigmaE', 'Range', 'Sigma',
                 'SigmaO', 'lamdba'),]
res.long <- melt(res, id.vars=c('knots', 'sd', 'par'))
g <- ggplot(res.long, aes(knots, value, group=par)) +
  facet_wrap('par', scales='free_y') + geom_line()
ggsave('plots/results_by_resolution.png', g, width=ggwidth, height=ggheight)
runtimes <- do.call(rbind, lapply(fits.res, function(x)
  data.frame(knots=x$n_knots, runtime=x$runtime)))
g <- ggplot(runtimes, aes(knots, runtime)) + geom_line()
ggsave('plots/runtime_by_resolution.png', g, width=ggwidth, height=ggheight)
}
### ------------------------------------------------------------
### OLD CODE

## ### Fit to some simulated data. Base if off the full model results from 3A.
## fit <- readRDS('results/..'
## d <- droplevels(subset(data.full, regcde=='3A'))
## fits.sim <- ldply(1:20, function(i)
##   simulate.fit(i=i, d=d, fit=fit, knots=1000, model='ST'))
## saveRDS(fits.sim, file='results/fits.sim.RDS')
## ggplot(fits.sim, aes(year, rel.error, group=rep)) + geom_line() + facet_wrap('trend')


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



