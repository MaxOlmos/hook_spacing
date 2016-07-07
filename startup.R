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


## Create mesh from real data
create.spde <- function(n_knots, make.plot=FALSE){
  ## df is used gobally here, so be careful
  loc_xy <- data.frame(x=df$longitude, y=df$latitude)
  knots <- kmeans( x=loc_xy, centers=n_knots )
  loc_centers <- knots$centers
  loc_xy <- cbind(loc_xy, cluster=knots$cluster, year=as.numeric(df$year))
  ## ggplot(loc_xy, aes(x,y, col=factor(cluster))) + geom_point(size=.5)
  mesh <- inla.mesh.create( loc_centers, refine=TRUE)
  spde <- inla.spde2.matern( mesh )
  if(make.plot){
    png(paste0('plots/mesh_', n_knots, '.png'), width=7, height=4, units='in', res=500)
    par(mar=.5*c(1,1,1,1))
    plot(mesh, main=NA, edge.color=gray(.7))
    points( df$longitude, df$latitude, cex=1, pch='.', col=rgb(0,0,0,.3))
    points( loc_centers, cex=.5, pch=20, col="red")
    dev.off()
  }
  return(list(mesh=mesh, spde=spde, cluster=knots$cluster))
}

make.inputs <- function(n_knots, model, likelihood=1, n_points_area=1e4, ...){
  spde <- create.spde(n_knots=n_knots)
  ## Calculate area of each SPDE grid by generating random points and
  ## seeing which proportion fall into each grid. This converges to area as
  ## the points go to Inf.
  ## loc_extrapolation <-
  ##   expand.grid(
  ##     "longitude"=seq(min(df$longitude), max(df$latitude),length=n_points_area),
  ##     "latitude"=seq(min(df$latitude), max(df$latitude),length=n_points_area))
  ## NN_extrapolation <- nn2( data=SimList$loc_xy, query=loc_extrapolation, k=1 )
  ## a_s <- table(factor(NN_extrapolation$nn.idx,levels=1:nrow(SimList$loc_xy))) / nrow(loc_extrapolation)
  ## Make inputs for all three models
  model <- match.arg(model, choices=c('NS', "S", "ST"))
  Data <- list(likelihood=likelihood, cph_i=df$cph,
               n_t=length(unique(df$year)),
               n_ft=max(df$spacing),
               s_i=spde$cluster-1,
               spacing_i=df$spacing,
               year_i=as.numeric(df$year)-1,
               month_i=as.numeric(df$month)-1,
               hooksize_i=as.numeric(df$hooksize)-1,
               geartype_i=as.numeric(df$geartype)-1,
               depth_i=df$depth,
               M0=spde$spde$param.inla$M0, M1=spde$spde$param.inla$M1,
               M2=spde$spde$param.inla$M2)
  Params <- list(intercept=5,
                 beta_year=rep(0, length(levels(df$year))),
                 beta_geartype= c(.17, .3, .3),
                 beta_month=rep(0, length(levels(df$month))),
                 beta_hooksize=rep(0, length(levels(df$hooksize))),
                 beta_depth=0, ln_tau_O=-.6, ln_tau_E=.25,
                 ln_kappa=.3,  ln_obs=-.2,
                 ln_spacing=0, spacing_devs=rep(0, length=Data$n_ft),
                 omega_s=rep(0,spde$mesh$n),
                 epsilon_st=matrix(0,nrow=nrow(Data$M0),ncol=Data$n_t))
  Random <- list(NS=c('spacing_devs'),
                 S=c('spacing_devs', 'omega_s'),
                 ST=c("omega_s", "epsilon_st", "spacing_devs"))
  ## Need to fix first level of each factor at 0 so they are
  ## identifiable. Get merged into the intercept. I.e., contrasts in R.
  list.factors <- list(
    beta_year=factor(c(NA, 1:(length(levels(df$year))-1))),
    beta_geartype=factor(c(NA, 1:(length(levels(df$geartype))-1))),
    beta_month=factor(c(NA, 1:(length(levels(df$month))-1))),
    beta_hooksize=factor(c(NA, 1:(length(levels(df$hooksize))-1))))
  Map <- list(
    ## Turn off spatial and spatio-temporal effects
    NS=c(list.factors, list(
      ln_tau_O = factor(NA), ln_kappa=factor(NA), ln_tau_E=factor(NA),
      omega_s=factor(rep(NA, len=spde$mesh$n)),
      epsilon_st=factor(rep(NA, len=spde$mesh$n*Data$n_t)))),
    ## Turn off just spatiotemporal effects
    S=c(list.factors, list(ln_tau_E=factor(NA),
      omega_s=factor(rep(NA, len=spde$mesh$n)),
      epsilon_st=factor(rep(NA, len=spde$mesh$n*Data$n_t)))),
    ## Turn on everything
    ST=list.factors)
      ## beta_statarea=factor(rep(NA, len=length(levels(df$statarea)))))))
  return(list(Data=Data, Params=Params, Random=Random[[model]], Map=Map[[model]]))
}

simulate.data <- function(loc_xy, loc_centers, params, n_years, SD_obs=.5,
                          Scale=1, SD_omega=1.5, SD_epsilon=.2 ){
    ## I'm passing locations into this instead of randomly generating them
    ## like Jim did.
    RF_omega <- RMgauss(var=SD_omega^2, scale=Scale)
    RF_epsilon <- RMgauss(var=SD_epsilon^2, scale=Scale)
    Omega_s <- RFsimulate(model=RF_omega, x=loc_centers[,1], y=loc_centers[,2])@data[,1]
    ## Generate spatio-temporal effects. Each year has a different number
    ## of observations so need to generate them separately and then
    ## carefully merge them in to add their effect
    Epsilon_st <- matrix(NA, nrow=nrow(loc_centers), ncol=n_years)
    for(t in 1:n_years){
      Epsilon_st[,t] <- RFsimulate(model=RF_epsilon,
                                   x=loc_centers[,1], y=loc_centers[,2])@data[,1]
    }
    x <- loc_xy
    x$epsilon_i <- Epsilon_st[cbind(x$cluster, x$year)]
    x$omega_i <- Omega_s[loc_xy$cluster]
    ## ggplot(x, aes(x,y, col=Omega_s)) + geom_point(size=.5, alpha=.5) + col.scale
    ## with(x, plot(cluster, Omega_s))
    ## y~1+year+geartype+month+hooksize+depth + spatial +spatiotemporal
    x$mu_i <- params$intercept + params$year[df$year]+
        params$geartype[df$geartype] +
        params$month[df$month]+params$hooksize[df$hooksize]+
        params$depth*df$depth + x$omega_i+
        x$epsilon_i
    ## Simulate samples for each site and year
    x$s_i <- x$cluster
    x$cph.simulated <- rnorm(n=length(x$mu_i), mean=x$mu_i, sd=SD_obs)
    return( x )
}