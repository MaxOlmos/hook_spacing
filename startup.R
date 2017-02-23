## Libraries used
library(TMB)
library(INLA)
library(RandomFields)
library(RANN)
library(ggplot2)
library(plyr)
library(reshape2)
library(maps)
library(maptools)
library(mapdata)
library(lubridate)
## data(wordHiresMapEnv)

## Global settings, mostly for plots
ggwidth <- 7
ggheight <- 4
mapdata <- map_data("world")
long.min <- -190; long.max <- -120
lat.min <- 40; lat.max <- 65
map.width <- 7*1.3
map.height <- 4*1.3
map.alpha <- .5
map.size <- .25
colors <- c("blue","green","red")
axis.col <- gray(.5)
effective.skates <- function(hooks.per.skate, hook.spacing, skates=1){
    (skates*hooks.per.skate)*1.52*(1-exp(-.06*hook.spacing))/100
}
effective.skates(100, 18)               # doesn't quite match up due to rounding
hook_power <- function(spacing, alpha, beta, lambda){
    (1- exp(-beta*spacing)^lambda ) / (1- exp(-beta*18)^lambda )
}

## This function will clean up the compiled files from a TMB object. Useful
## for development when need to ensure all change are propogated.
clean.TMB.files <- function(model.path){
  o <- paste0(model.path,'.o')
  dll <- paste0(model.path, '.dll')
  tryCatch(expr=dyn.unload(dynlib(model.path)), error=function(cond){x=3})
  if(file.exists(dll)) trash <- file.remove(dll)
  if(file.exists(o)) trash <- file.remove(o)
}
## Forces elements of x to be unique by adding number to any duplicated
## entries. Helper function for plotting functions below
add.unique.names <- function(x){
  as.vector(unlist(llply(unique(x), function(y){
  z <- x[x %in% y]
  if(length(z)>1) z <- paste0(z,"_", seq_along(z))
  z})))
}

## Pass it a list of fits from run.logbook, and plots parameter comparisons
## between the two
plot.parameter.comparison <- function(fits, level.name, levels){
  test <- ldply(1:length(fits), function(x)
    cbind(level.name=levels[x], fits[[x]]$sd.par))
  test[,level.name] <- test$level.name
  test$level.name <- NULL
  test$par3 <- as.numeric(as.factor(paste0(test$par2, "_", test[,level.name])))
  g <- ggplot(test, aes_string('par3', 'value', color=level.name,
  group=level.name, ymin='lwr', ymax='upr')) +
    geom_linerange(lwd=1.5) + facet_wrap('par', scales='free') +
    geom_hline(yintercept=0) + geom_point()
  return(g)
}

## Same as above but for spacing effect. Currently broken.
plot.power.comparison <-
  function(fits){
    test <- ldply(1:length(fits), function(x)
      cbind(model=fits[[x]]$model.name, form=fits[[x]]$form.name, fits[[x]]$sd.hook_power))
    g <- ggplot(test, aes(hook_power, value, fill=model, col=model,
                          ymin=lwr, ymax=upr)) +
      geom_ribbon(lwd=1, alpha=1/3) + facet_grid(form~.) +
      geom_vline(xintercept=18) + geom_line(lwd=1) + ylab("Hook Power")
    return(g)
  }

## Same as above but for abundance trends
plot.cpue.comparison <- function(fits){
  test <- ldply(1:length(fits), function(x)
    cbind(model=fits[[x]]$model.name, form=fits[[x]]$form.name, year=1996:2015, fits[[x]]$sd.density))
  g <- ggplot(test, aes(year, value, color=model, fill=model, group=model, ymin=lwr, ymax=upr)) +
    geom_ribbon(alpha=1/3) + facet_grid(form~., scales='free_y') + geom_line(lwd=2)
  g <- g + theme_bw() + ylab("Relative Abundance")
  return(g)
}

## Same as above but for residuals
plot.resids.comparison <- function(fits){
  test <- ldply(1:length(fits), function(x)
    data.frame(model=fits[[x]]$model.name, form=fits[[x]]$form.name,
               resids=fits[[x]]$report$resids, preds=fits[[x]]$report$preds))
  g <- ggplot(test, aes(x=resids, group=model, fill=model)) +
    facet_grid(form~.) + geom_histogram(alpha=1/3, bins=200, position='identity')
  g <- g + theme_bw()
  return(g)
}


## Create mesh from real data
create.spde <- function(dat, n_knots, make.plot=FALSE, jitter=.3){
  loc_xy <- data.frame(x=dat$longitude, y=dat$latitude)
  knots <- kmeans( x=loc_xy, centers=n_knots )
  loc_centers <- knots$centers
  loc_xy <- cbind(loc_xy, cluster=knots$cluster, year=as.numeric(dat$year))
  ## ggplot(loc_xy, aes(x,y, col=factor(cluster))) + geom_point(size=.5)
  mesh <- inla.mesh.create( loc_centers, refine=TRUE)
  spde <- inla.spde2.matern( mesh )
  if(make.plot){
    png(paste0('plots/mesh_', n_knots, '.png'), width=7, height=4, units='in', res=500)
    par(mar=.5*c(1,1,1,1))
    plot(mesh, main=NA, edge.color=gray(.7))
    points( jitter(dat$longitude, amount=jitter), jitter(dat$latitude, amount=jitter), cex=1, pch='.', col=rgb(0,0,0,.3))
    points( loc_centers, cex=.5, pch=20, col="red")
    dev.off()
  }
  return(list(mesh=mesh, spde=spde, cluster=knots$cluster, loc_centers=loc_centers))
}

make.inputs <- function(dat, n_knots, model, form,
  likelihood=1, vessel_effect=FALSE,  n_points_area=1e3, ...){
  spde <- create.spde(dat=dat, n_knots=n_knots)
  loc <- spde$mesh$loc[,1:2]
  n_s <- nrow(loc)
  ## Calculate area of each SPDE grid by generating random points and
  ## seeing which proportion fall into each grid. This converges to area as
  ## the points go to Inf.
  loc_extrapolation <-
    expand.grid(
      "longitude"=seq(min(dat$longitude), max(dat$longitude),length=n_points_area),
      "latitude"=seq(min(dat$latitude), max(dat$latitude),length=n_points_area))
  NN_extrapolation <- nn2( data=loc, query=loc_extrapolation, k=1 )
  area_s <- table(factor(NN_extrapolation$nn.idx,levels=1:n_s)) / nrow(loc_extrapolation)
  ## hist(area_s)
  ## test <- cbind(loc_extrapolation, idx=NN_extrapolation$nn.idx)
  ## ggplot(test, aes(longitude,latitude, color=idx)) + geom_point()

  ## Make inputs for all three models
  model <- match.arg(model, choices=c('NS', "S", "ST"))
  if(model=='NS') space <- 0
  if(model=='S') space <- 1
  if(model=='ST') space <- 2
  Data <- list(likelihood=likelihood, form=form, space=space,
               catch_i=dat$catch,
               hooks_i=dat$hooks,
               n_t=length(unique(dat$year)),
               n_s=n_s,
               n_ft=max(dat$spacing)+5,
               n_v=length(unique(dat$vessel)), # no. vessels
               s_i=spde$cluster-1,
               spacing_i=dat$spacing,
               year_i=as.numeric(dat$year)-1,
               month_i=as.numeric(dat$month)-1,
               hooksize_i=as.numeric(dat$hooksize)-1,
               geartype_i=as.numeric(dat$geartype)-1,
               vessel_i=as.numeric(dat$vessel)-1,
               depth_i=dat$depth,
               area_s=area_s,
               M0=spde$spde$param.inla$M0, M1=spde$spde$param.inla$M1,
               M2=spde$spde$param.inla$M2)
  Params <- list(intercept=.5,
                 beta_year=rep(0, length(levels(dat$year))),
                 beta_geartype= c(0,0, 0),
                 beta_month=rep(0, length(levels(dat$month))),
                 beta_hooksize=rep(0, length(levels(dat$hooksize))),
                 beta_depth=0,beta_depth2=0,
                 beta_spacing=0.5, lambda=1,
                 ln_tau_O=-.6, ln_tau_E=.25,
                 ln_kappa=.3,  ln_obs=-.2, ln_spacing=0,
                 ln_vessel=.1,
                 spacing_devs=rep(0, length=Data$n_ft),
                 vessel_v=rep(0, length=length(unique(Data$vessel))),
                 omega_s=rep(0,spde$mesh$n),
                 epsilon_st=matrix(0,nrow=nrow(Data$M0),ncol=Data$n_t))

  ## If vessel effect turned on, set it to random
  v <- if(vessel_effect) 'vessel_v' else NULL
  if(form==1){
    Random <- list(NS=c('spacing_devs', v),
                   S=c('spacing_devs', 'omega_s', v),
                   ST=c("omega_s", "epsilon_st", "spacing_devs",v))
  } else { ## turn off spacing RE if using H&S formula
    Random <- list(NS=v,
                   S=c('omega_s',v),
                   ST=c("omega_s", "epsilon_st",v))
  }
  ## Need to fix first level of each factor at 0 so they are
  ## identifiable. Get merged into the intercept. I.e., contrasts in R.
   list.factors <- list(
     beta_year=factor(c(NA, 1:(length(levels(dat$year))-1))),
    beta_geartype=factor(c(NA, 1:(length(levels(dat$geartype))-1))),
    ## beta_month=factor(c(NA, 1:(length(levels(dat$month))-1))),
    ## beta_hooksize=factor(c(NA, 1:(length(levels(dat$hooksize))-1))),
    beta_month=factor(rep(NA, length(levels(dat$month)))),
    beta_hooksize=factor(rep(NA, length(levels(dat$hooksize)))))

  ## Turn off parameters for spacing depending on the form
  if(form==1) {
    ## random walk on spacing
    xx <- list(beta_spacing=factor(NA), lambda=factor(NA))
  }
  if(form==2) {
    ## H&S form on spacing
    xx <- list(spacing_devs=factor(rep(NA, length=Data$n_ft)),
               ln_spacing=factor(NA))
  }
  if(form==3) {
    ## No effect (set to zero in the TMB model)
    xx <- list(spacing_devs=factor(rep(NA, length=Data$n_ft)),
               ln_spacing=factor(NA), beta_spacing=factor(NA),
                lambda=factor(NA))
  }
  ## add vessels if needed
  if(!vessel_effect) xx <- c(xx, list(vessel_v=factor(rep(NA, length=Data$n_v))), list(ln_vessel=factor(NA)))
  Map <- list(
    ## Turn off spatial and spatio-temporal effects
    NS=c(list.factors, xx, list(
      ln_tau_O = factor(NA), ln_kappa=factor(NA), ln_tau_E=factor(NA),
      omega_s=factor(rep(NA, len=spde$mesh$n)),
      epsilon_st=factor(rep(NA, len=spde$mesh$n*Data$n_t)))),
    ## Turn off just spatiotemporal effects
    S=c(list.factors, xx, list(ln_tau_E=factor(NA),
      epsilon_st=factor(rep(NA, len=spde$mesh$n*Data$n_t)))),
    ## Turn on everything
    ST=c(xx,list.factors))
  return(list(Data=Data, Params=Params, Random=Random[[model]], Map=Map[[model]]))
}

model.name <- function(model) switch(model, NS="No Space", S='Space',
                       ST='Spatiotemporal')
form.name <- function(form)  c('Smoother', 'Parametric', 'Constant')[form]


run.logbook <- function(dat, n_knots, model, form, vessel_effect, likelihood=1, trace=10){

  message(paste0('Starting model ', model, ', form=',form, ', vessel=', vessel))
  ## Causes a nasty bug if extra factor levels are around so drop those
  dat <- droplevels(dat)
  n_years <- length(unique(dat$year))
  model.name <- model.name(model)
  form.name <- form.name(form)
  start <- Sys.time()
  Inputs <- make.inputs(dat=dat, n_knots=n_knots, model=model, form=form,
                        vessel_effect=vessel_effect, likelihood=likelihood)
  Obj <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params,
                   random=Inputs$Random, map=Inputs$Map)
  Obj$env$beSilent()
  ## Set bounds for parameters via limits in optimizer
  lower <- rep(-Inf, len=length(unlist(Inputs$Params)))
  upper <- rep(Inf,  len=length(unlist(Inputs$Params)))
  names(upper) <- names(lower) <- names(unlist(Inputs$Params))
  lower['gamma'] <- 0; upper['gamma'] <- 5
  Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
                lower=lower, upper=upper,
                control=list(trace=trace, eval.max=1e4, iter.max=1e4 ))
  report.temp <- Obj$report();
  sd.temp <- sdreport(Obj)
  sd.df <- data.frame(par=names(sd.temp$value), value=sd.temp$value, sd=sd.temp$sd)
  sd.df <- within(sd.df,{
    upr <- value+1.96*sd
    lwr <- value-1.96*sd
    table <- paste0(format(value,digits=3, scientific=TRUE),
                  ' (',format(lwr,digits=3, sci=TRUE), '-',
                  format(upr,digits=3, sci=TRUE),')')
  })
  sd.df$par <- as.character(sd.df$par)
  sd.df$par2 <- add.unique.names(sd.df$par)
  ## sd.spacing <- sd.df[sd.df$par=='spacing',]
  ## sd.spacing$spacing <- seq_along(sd.spacing$par)
  sd.hook_power <- sd.df[grep('hook_power', x=sd.df$par),]
  sd.hook_power$hook_power <- seq_along(sd.hook_power$par)
  sd.cpue <- sd.df[grep('cph_t', x=sd.df$par),]
  sd.density <- sd.df[grep('area_weighted_density_t', x=sd.df$par),]
  sd.density$par <- paste0('density_', 1:n_years)
  sd.cpue$year <- sd.density$year <- 1996:2015
  sd.par <- sd.df[-grep('hook_power|cph_t|area_weighted_density_t', x=sd.df$par),]
  ## This is the order of pars I want to output for easy conversion into a table.
  table.pars <-
    c('intercept', 'beta_depth', 'beta_depth2', 'beta_geartype', 'Sigma',
                'Sigma_vessel', 'Range', 'SigmaE', 'SigmaO')
  sd.table <- subset(sd.par, par %in% table.pars, select=c(par, table))
  runtime <- as.numeric(difftime(Sys.time(),start, units='mins'))
  x <- list(model=model, n_knots=n_knots, form=form, likelihood=likelihood,
            model.name=model.name, form.name=form.name,
            runtime=runtime, report=report.temp,
                sd.hook_power=sd.hook_power, sd.table=sd.table,
            sd.density=sd.density, sd.cpue=sd.cpue, sd.par=sd.par, Obj=Obj, Opt=Opt,
            Inits=Inputs$Params, Map=Inputs$Map, Random=Inputs$Random)
  return(x)
}

#' @param data A data frame with the original data
#' @param fit A fit to  an existing model. The parameter estimates will be
#'   used.
#' @param beta A coefficient for altering the spacing. If NULL the original
#'   spacing left.
#' @return A modified data set where catch has been replaced
simulate.data <- function(dat, n_knots, fit, beta=NULL){
  dat <- droplevels(dat)
  ## Create simulated spatial process, using output from fit
  n_points_area <- 1e3
  n_years <- length(unique(dat$year))
  spde <- create.spde(dat=dat, n_knots=n_knots)
  loc <- spde$mesh$loc[,1:2]
  loc_centers <- spde$loc_centers
  n_s <- nrow(loc)
  ## Calculate area of each SPDE grid by generating random points and
  ## seeing which proportion fall into each grid. This converges to area as
  ## the points go to Inf.
  loc_extrapolation <-
    expand.grid(
      "longitude"=seq(min(dat$longitude), max(dat$longitude),length=n_points_area),
      "latitude"=seq(min(dat$latitude), max(dat$latitude),length=n_points_area))
  NN_extrapolation <- nn2( data=loc, query=loc_extrapolation, k=1 )
  area_s <- table(factor(NN_extrapolation$nn.idx,levels=1:n_s)) / nrow(loc_extrapolation)
  SD_omega <- fit$sd.par$value[fit$sd.par$par=='SigmaO']
  SD_epsilon <- fit$sd.par$value[fit$sd.par$par=='SigmaE']
  Sigma <- fit$sd.par$value[fit$sd.par$par=='Sigma']
  Scale <- .2#fit$sd.par$value[fit$sd.par$par=='Range']
  intercept <- fit$sd.par$value[fit$sd.par$par=='intercept']
  beta_years <- fit$sd.par$value[grep('beta_year', fit$sd.par$par2)]
  beta_geartype <- fit$sd.par$value[grep('beta_geartype', fit$sd.par$par2)]
  beta_spacing <- fit$sd.par$value[fit$sd.par$par=='beta_spacing']
  lambda <- fit$sd.par$value[fit$sd.par$par=='lambda']

  ## Spatial effects
  RF_omega <- RMgauss(var=SD_omega^2, scale=Scale)
  Omega_s <- RFsimulate(model=RF_omega, x=loc_centers[,1], y=loc_centers[,2])@data[,1]

  ## Spatio-temporal effects. Each year has a different number
  ## of observations so need to generate them separately and then
  ## carefully merge them in to add their effect
  RF_epsilon <- RMgauss(var=SD_epsilon^2, scale=Scale)
  log_D_st <- Epsilon_st <- matrix(NA, nrow=nrow(loc_centers), ncol=n_years)
  for(t in 1:n_years){
    Epsilon_st[,t] <-
      RFsimulate(model=RF_epsilon, x=loc_centers[,1], y=loc_centers[,2])@data[,1]
    log_D_st[,t] <- intercept + Epsilon_st[,t] + Omega_s + beta_years[t]
  }
  density_t = colSums( exp(log_D_st) )
  ## Now match up expected density with each set
  dat$density <- exp(log_D_st[cbind(spde$cluster, as.numeric(dat$year))])

  ## ggplot(dat, aes(longitude,latitude, col=cluster)) +
  ##   geom_point(size=.5, alpha=.5) + scale_color_gradient(low='blue', high='red')
  ## Now sample from the truth
  ## y~1+year+geartype+month+hooksize+depth + spatial +spatiotemporal

  ## In some cases want to alter the spacing to create trends in hook
  ## spacing. Assuming for now a simple linear decrease over time. If
  ## beta=1 this declines, if beta=0 it is constant. If NULL the original
  ## spacing is used.
  if(!is.null(beta)){
    years <- as.numeric(dat$year)
    dat$spacing <- pmax(3, rnbinom(n=nrow(dat), size=35,
                                    mu=20-beta*years))
    ##    dat$spacing <- pmax(3, rpois(n=nrow(dat), lambda=15-beta*years))
     ggplot(dat, aes(years, spacing)) + geom_jitter(alpha=.05) + geom_smooth()
  }
  ## Expected catch is density*q*hooks*f(spacing)
  dat$mu_i <-
    dat$density*                       # denisty
    exp(beta_geartype[dat$geartype])*  # q
    dat$hooks *                        # hooks
    hook_power(dat$spacing, alpha=1, beta=beta_spacing, lambda=lambda) # hook power

  ## Simulate samples for each site and year
  dat$catch <- exp(rnorm(n=length(dat$mu_i), mean=log(dat$mu_i), sd=Sigma))
  return( list(data=dat, density_t=density_t) )
}

plot.spacing.fit <- function(results, multiple.fits=FALSE){
  if(multiple.fits){
    results <- ldply(results, function(y) cbind(model=y$model, likelihood=y$likelihood, form=y$form,
                                                y$sd.spacing))
  } else {
    results <- results$sd.spacing
  }
  print(str(results))
  g <- ggplot(results, aes(spacing, value, ymax=value+2*sd, ymin=value-2*sd)) +
    geom_pointrange()
  if(multiple.fits) g <-  g+facet_grid(model+likelihood~form)
  g
}

make.inputs.experimental <- function(hs){
  Data <-
    list(n_i=nrow(hs), n_s=length(unique(hs$group)),
         catch_i=hs$catch, hooks_i=hs$hooks, day_i=hs$daynumber,
         spacing_i=hs$spacing, site_i=hs$group-1)
  Params <-
    list(eta_mean=.3, eta_sd=.1, sigma_mean=.5, sigma_sd=.2,
         gamma=0.05, beta=.05, lambda=1,
         eta_s=rep(.3, len=Data$n_s), sigma_s=rep(.5,len=Data$n_s))
  Map <- list(lambda=factor(NA))          # can't estimate this so fix at 1
  ## Set bounds for parameters via limits in optimizer
  lower <- rep(-Inf, len=length(unlist(Params)))
  upper <- rep(Inf,  len=length(unlist(Params)))
  names(upper) <- names(lower) <- names(unlist(Params))
  lower['gamma'] <- 0; upper['gamma'] <- 1
  lower['lambda'] <- 0; upper['lambda'] <- 5
  lower['eta_sd'] <- 0; upper['eta_sd'] <- Inf
  lower['sigma_sd'] <- 0; upper['sigma_sd'] <- Inf
  Random <- c('eta_s', 'sigma_s')
  Inputs <- list(Data=Data, Params=Params, Random=Random, Map=Map, Lower=lower, Upper=upper)
}


simulate.fit <- function(i, d, fit, knots, model){
  ## ## First simulate with original spacing levels
  ## sim.orig <- simulate.data(d, knots, fit=fit, beta=NULL)
  ## fit.orig <- run.logbook(dat=sim.orig$data, n_knots=knots, model=model, form=3, vessel=FALSE)
  ## Now resimulate with a totally flat hook spacing
  print(i)
  sim.flat <- simulate.data(d, knots, fit=fit, beta=0)
  fit.flat <- run.logbook(sim.flat$data, n_knots=knots, model=model,
                          form=3, vessel=FALSE, trace=100)
  temp.flat <- data.frame(rep=i, trend='Flat', year=1996:2015, density=sim.flat$density_t)
  ## And again with a decreasing trend in spacing over time
  sim.trend <- simulate.data(d, knots, fit=fit, beta=.5)
  fit.trend <- run.logbook(sim.trend$data, n_knots=knots, model=model,
                           form=3, vessel=FALSE, trace=100)
  temp.trend <- data.frame(rep=i, trend='Trend', year=1996:2015, density=sim.trend$density_t)
  ## Return only what needed
  res <- rbind(
    ## data.frame(rep=i, trend='Original', year=1996:2015,
    ##            rel.abundance.true=sim.orig$density_t/mean(sim.orig$density_t),
    ##            rel.abundance.fit=fit.orig$sd.density$value/mean(fit.orig$sd.density$value)),
    data.frame(rep=i, trend='Flat', year=1996:2015,
               rel.abundance.true=sim.flat$density_t/mean(sim.flat$density_t),
               rel.abundance.fit=fit.flat$sd.density$value/mean(fit.flat$sd.density$value)),
    data.frame(rep=i, trend='Trend', year=1996:2015,
               rel.abundance.true=sim.trend$density_t/mean(sim.trend$density_t),
               rel.abundance.fit=fit.trend$sd.density$value/mean(fit.trend$sd.density$value)))
  if(i==1) saveRDS(list(sim.flat=sim.flat, sim.trend=sim.trend), file='results/data.sim.RDS')
  res <- within(res, {
    abs.error <- -rel.abundance.true+rel.abundance.fit
    rel.error <- (-rel.abundance.true+rel.abundance.fit)/rel.abundance.true})
  return(res)
}
