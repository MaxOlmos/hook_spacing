
### ------------------------------------------------------------
## Step 2. Run models. Models= no spatial effect (NS), spatial model (S)
## and full spatio-temporal (ST)
Version <- "spatiotemporal_cpue_spacing"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
## ## Test model is working
## Inputs <- make.inputs(n_knots=20, model='ST', likelihood=1)
## Obj <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params,
##                  random=Inputs$Random, map=Inputs$Map)
## Obj$fn()

## Loop through each model, running and saving results.
n_knots <- 500
Opt.list <- Report.list <- SD.list <- list()
for(m in c('NS', "S", "ST")){
  print(paste0('Starting model: ', m))
  Inputs <- make.inputs(n_knots=n_knots, model=m, likelihood=1)
  Obj <- MakeADFun( data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random,
                   map=Inputs$Map)
  trash <- Obj$env$beSilent()
  start <- Sys.time()
  temp <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
                 control=list(trace=50, eval.max=1e4, iter.max=1e4))
  Opt.list[[m]] <- nlminb( start=temp$par, objective=Obj$fn, gradient=Obj$gr,
                     control=list(trace=50, eval.max=1e4, iter.max=1e4))
  Opt.list[[m]][["final_diagnostics"]] <-
    data.frame( "Name"=names(Obj$par), "final_gradient"=as.numeric(Obj$gr(Opt.list[[m]]$par)))
  Report.list[[m]] <- test <- Obj$report()
  Report.list[[m]]$time <- as.numeric(difftime(Sys.time(),start,
                                               units='mins'))
  SD.list[[m]] <- sdreport(Obj)
  ##  make.model.plots(

}
saveRDS(Report.list, 'report_models.RDS')
saveRDS(SD.list, 'sd_models.RDS')

## Process the results and make plots
Report.list <- readRDS('report_models.RDS')
Results <- do.call(rbind, lapply(Report.list, function(x)
  data.frame(x[c('intercept', 'beta_depth','SigmaE', 'Range', 'Sigma',
                 'SigmaO', 'jnll', 'time')])))
Results$model <- c('M1', "M2", "M3")
SD.list <- readRDS('sd_models.RDS')
cpue.trends <- (lapply(SD.list, function(x)
  data.frame(year=1995:2012, cpue=x['value'], se=x['sd'])))
cpue.trends <- melt(cpue.trends, id.vars=c('year', 'value', 'sd'))
names(cpue.trends) <- c('year', 'cpue', 'sd', 'model')
cpue.trends <- within(cpue.trends, {
                        lwr=cpue-sd
                        upr=cpue+sd})

## CPUE trends with and without errors
g <- ggplot(cpue.trends, aes(year, cpue, group=model, color=model)) +
  geom_line() + ylab("Relative CPUE")
ggsave('plots/cpue_trends.png', g, width=ggwidth, height=ggheight)
g <- ggplot(cpue.trends, aes(year, cpue, ymin=lwr, ymax=upr)) +
  geom_linerange(color=gray(.5)) + geom_point()+ylab("Relative CPUE") + facet_wrap('model')
ggsave('plots/cpue_trends_errorbars.png', g, width=ggwidth, height=ggheight)

## Model parameters and meta data
Results.long <- melt(Results, id.vars='model')
g <- ggplot(Results.long, aes(model, (value), group=variable)) +
  facet_wrap('variable', scales='free_y') + geom_line()
ggsave('plots/results_by_model.png', g, width=ggwidth, height=ggheight)
qqplots <- ldply(Report.list, function(x){
                   xx <- qqnorm(y=x$resids, plot=FALSE)
                 return(data.frame(xx))})
g <- ggplot(qqplots, aes(x,y)) + geom_point() + facet_wrap('.id') + geom_abline(slope=1,intercept=0)
ggsave('plots/qqresids.png', g, width=ggwidth, height=ggheight)

## Plots of resids to see if there's a spatial trend removed
df$resids.M1 <- Report.list[['M1']]$resids
df$resids.M2 <- Report.list[['M2']]$resids
df$resids.M3 <- Report.list[['M3']]$resids
col.scale <-  scale_colour_distiller(palette='Spectral', limits=c(-3,3)) #scale_color_continuous(low='blue', high='red')
myxlim <- xlim(-155,-145); myylim <- ylim(56,61);
## Some global ones
g <- ggplot(df, aes(longitude, latitude, color=resids.M1)) +
  geom_point(alpha=1, size=.15) + col.scale
ggsave('plots/spatial_errors_M1.png',g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(longitude, latitude, color=resids.M2)) +
  geom_point(alpha=1, size=.15) + col.scale
ggsave('plots/spatial_errors_M2.png',g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(longitude, latitude, color=resids.M3)) +
  geom_point(alpha=1, size=.15) + col.scale
ggsave('plots/spatial_errors_M3.png',g, width=ggwidth, height=ggheight)

for(yr in unique(df$year)){
g <- ggplot(subset(df, year==yr), aes(longitude, latitude, color=resids.M1)) +
  geom_point(alpha=1, size=.5) + col.scale + myxlim+myylim
ggsave(paste0('plots/annual/spatial_errors_',yr,'_M1.png'),g, width=7, height=4)
g <- ggplot(subset(df, year==yr), aes(longitude, latitude, color=resids.M2)) +
  geom_point(alpha=1, size=.5) + col.scale + myxlim+myylim
ggsave(paste0('plots/annual/spatial_errors_',yr,'_M2.png'),g, width=7, height=4)
g <- ggplot(subset(df, year==yr), aes(longitude, latitude, color=resids.M3)) +
  geom_point(alpha=1, size=.5) + col.scale + myxlim+myylim
ggsave(paste0('plots/annual/spatial_errors_',yr,'_M3.png'),g, width=7, height=4)
}


### OLD CODE
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

## Obj$env$beSilent()
## temp <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
##                  control=list(trace=10, eval.max=1e4, iter.max=1e4))
## report.temp <- Obj$report()
## temp2 <- data.frame(resids=report.temp$resids, spacing=df$spacing)
## plot(report.temp$spacing_std)
## plot(report.temp$spacing)
## hist(report.temp$spacing)
## sd.temp <- sdreport(Obj)
## value.spacing <- sd.temp$value[grep('spacing_std', x=names(sd.temp$value))]
## sd.spacing <- sd.temp$sd[grep('spacing_std', x=names(sd.temp$value))]
## df.spacing <- data.frame(ft=seq_along(value.spacing), value=value.spacing, sd=sd.spacing)[-1,]
## ggplot(df.spacing, aes(ft, value, ymin=value-2*sd, ymax=value+2*sd)) + geom_errorbar()

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
