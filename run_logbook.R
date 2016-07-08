### ------------------------------------------------------------
## Step 2. Run models. Models= no spatial effect (NS), spatial model (S)
## and full spatio-temporal (ST)
Version <- "models/spatiotemporal_cpue_spacing"
dyn.unload( dynlib(Version) )
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

## Run ST model with and without the HS formula
for(form in 1:2){
  Inputs <- make.inputs(n_knots=1000, model='ST', form=form, likelihood=1)
  Obj <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params,
                   random=Inputs$Random, map=Inputs$Map)
  Obj$env$beSilent()
  Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
                control=list(trace=10, eval.max=1e4, iter.max=1e4))
  report.temp <- Obj$report(); sd.temp <- sdreport(Obj)
  value.spacing <- sd.temp$value[grep('spacing_std', x=names(sd.temp$value))][-1]
  sd.spacing <- sd.temp$sd[grep('spacing_std', x=names(sd.temp$value))][-1]
  uncertainty.df <- data.frame(spacing=seq_along(value.spacing), value=value.spacing, sd=sd.spacing)
  logbook.results <- list(Obj=Obj, Opt=Opt, report=report.temp,
                             uncertainty.df=uncertainty.df)
  if(form==1) saveRDS(logbook.results, file='results/logbook.results.RDS')
  if(form==2) saveRDS(logbook.results, file='results/logbook.hs.results.RDS')
  rm(sd.temp, report.temp, value.spacing, sd.spacing, uncertainty.df,
     logbook.results)
  gc(); gc()
}


### ------------------------------------------------------------
### OLD CODE
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
