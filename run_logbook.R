### ------------------------------------------------------------
## Step 1: read in and prep the data for this model
source('startup.R')
df.unfiltered <- readRDS(file='data/data_unfiltered.RDS')
df.unfiltered$spacing <- round(df.unfiltered$spacing)
df <- readRDS(file='data/data.RDS')
df$spacing <- round(df$spacing)
## df.simulated <- df
## df.simulated$catch <- with(df, hooks*exp(-.1+rnorm(n=nrow(df), mean=0, sd=.1)))
## df <- df.simulated
n_years <- length(unique(df$year))
Version <- "models/spatiotemporal_cpue_spacing"
clean.TMB.files(Version)
compile( paste0(Version,".cpp"))
dyn.load( dynlib(Version))

### ------------------------------------------------------------
## Step 2. Model exploration. Models= no spatial effect (NS), spatiaggl model (S)
## and full spatio-temporal (ST). orm=1 implies a random walk on hook
## spacing, form=2 is the parametric HS model.


### Explore effects of key dimensions.

## Spatial effect
knots <- 1000
ns <- run.logbook(n_knots=knots, model='NS', form=2, vessel=TRUE)
s <- run.logbook(n_knots=knots, model='S', form=2, vessel=TRUE)
st <- run.logbook(n_knots=knots, model='ST', form=2, vessel=TRUE)
g <- plot.parameter.comparison(list(ns,s,st),
     level.name='model', levels=c('NS', 'S', 'ST'))
ggsave('plots/par_comparison_model.png')

## Spacing effect for NS model
knots <- 10
ns1 <- run.logbook(n_knots=knots, model='NS', form=1, vessel=FALSE)
ns2 <- run.logbook(n_knots=knots, model='NS', form=2, vessel=FALSE)
ns3 <- run.logbook(n_knots=knots, model='NS', form=3, vessel=FALSE)
g <- plot.parameter.comparison(list(ns1, ns2, ns3),
     level.name='form', levels=c('Smoother', 'H&S', 'None'))
ggsave('plots/par_comparison_model.png', g, width=12, height=6)
g <- plot.spacing.comparison(fits=list(ns1, ns2, ns3))
ggsave('plots/spacing_comparison.png', g, width=5, height=6)
g <- plot.power.comparison(fits=list(ns1, ns2, ns3))
ggsave('plots/power_comparison.png', g, width=5, height=6)
## plot(df$spacing, ns2$report$resids, pch='.', col=rgb(0,0,0,.1))

## Spacing vs model!
knots <- 1000
form1 <- run.logbook(n_knots=knots, model='NS', form=1, vessel=TRUE)
form2 <- run.logbook(n_knots=knots, model='NS', form=2, vessel=TRUE)
form3 <- run.logbook(n_knots=knots, model='NS', form=3, vessel=TRUE)
form4 <- run.logbook(n_knots=knots, model='S', form=1, vessel=TRUE)
form5 <- run.logbook(n_knots=knots, model='S', form=2, vessel=TRUE)
form6 <- run.logbook(n_knots=knots, model='S', form=3, vessel=TRUE)
form7 <- run.logbook(n_knots=knots, model='ST', form=1, vessel=TRUE)
form8 <- run.logbook(n_knots=knots, model='ST', form=2, vessel=TRUE)
form9 <- run.logbook(n_knots=knots, model='ST', form=3, vessel=TRUE)
fits <- list(form1, form2, form3, form4, form5, form6, form7, form8, form9)
saveRDS(fits, file='results/fits_form_vs_model.RDS')
## Make quick plots
g <- plot.parameter.comparison(fits=list(form7,form8,form9),
 level.name='model', levels=c('Nonparametric', 'Hamley & Skud', 'None'))
ggsave('plots/par_comparison_form.png', g, width=10, height=6)
g <- plot.cpue.comparison(fits=fits)
ggsave('plots/cpue_comparison_form.png', g, width=10, height=6)
g <- plot.spacing.comparison(fits=fits)
ggsave('plots/spacing_comparison.png', g, width=5, height=6)
g <- plot.resids.comparison(fits=fits)
ggsave('plots/resids_comparison.png', g, width=5, height=6)

## Cleanup
dyn.unload( dynlib(Version))

### ------------------------------------------------------------
### OLD CODE

## ## Test code. Track parameter traces for plotting
## Inputs <- make.inputs(n_knots=50, model='S', form=1, likelihood=1)
## Obj <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params,
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
##   Obj <- MakeADFun( data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random,
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
##   g+geom_point(data=df, aes(x=longitude, y=latitude, colour=statarea),
##                size=map.size, alpha=map.alpha)
## ggsave('plots/map_statareas.png', g2, width=map.width, height=map.height)
## png('plots/pairs.png', width=ggwidth, height=ggheight, units='in', res=500)
## pairs(df[, c('logcpue', 'hooksize', 'geartype', 'depth', 'month')])
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


