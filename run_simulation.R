## Fit a quick linear model to get crude estimates of parameters for
## simulating data sets.
mod.lm <- lm(logcpue~year+geartype+month+hooksize+depth, data=df)
summary(mod.lm)
coefs <- coef(mod.lm)
x <- names(coefs)
params.sim <- list(
  year=c(0,as.numeric(coefs[grep('year', x=x)])),
  month=c(0,as.numeric(coefs[grep('month', x=x)])),
  hooksize=c(0,as.numeric(coefs[grep('hooksize', x=x)])),
  geartype=c(0,as.numeric(coefs[grep('geartype', x=x)])),
  intercept=as.numeric(coefs['(Intercept)']),
  depth=as.numeric(coefs['depth']))
params.sim <- list(
  year=rep(0, 18),
  month=rep(0,10),
  hooksize=rep(0,12),
  geartype=rep(0,3),
  intercept=as.numeric(coefs['(Intercept)']),
  depth=0)
rm(coefs, x, mod.lm)



params.sim <- list(
  year=rep(0, 18),
  month=rep(0,10),
  hooksize=rep(0,12),
  geartype=rep(0,3),
  intercept=as.numeric(coefs['(Intercept)']),
  depth=0)

logmean=1; Scale=0.2; SD_omega=1; SD_epsilon=1; n_per_year=100; n_years=18
simulate.data <- function(loc_xy, loc_centers, params, n_years, SD_obs=.00008,
                          Scale=5, SD_omega=1.5, SD_epsilon=.2 ){
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
    x$logcpue.simulated <- rnorm(n=length(x$mu_i), mean=x$mu_i, sd=SD_obs)
    return( x )
}

## Generate and plot the global spatial random effectsd
df.simulated <-
  simulate.data(loc_xy, loc_centers, params.sim, Scale=.5, SD_obs=.0001, SD_epsilon=.2,
                n_years=18)
col.scale <-  scale_colour_distiller(palette='RdYlGn') #scale_color_continuous(low='blue', high='red')
## g <- ggplot(subset(df.simulated, s_i<10), aes(x, y, color=factor(s_i))) +
##   geom_point(alpha=.5, size=.25)
## Add them on to master data.frame
df$logcpue.simulated <- df.simulated$logcpue.simulated
df$mu_i <- df.simulated$mu_i
df$omega_i <- df.simulated$omega_i
g <- ggplot(subset(df, year==1995), aes(longitude, latitude, color=omega_i)) +
  geom_point(alpha=.5, size=.25) + col.scale
ggsave('plots/spatial_errors.png',g, width=7, height=4)
g <- ggplot(subset(df, year==1995), aes(longitude, latitude, color=mu_i)) +
  geom_point(alpha=.5, size=.25) +
   col.scale #+ facet_wrap('year')
##      xlim(-155, -145) + ylim(58,60) + theme_bw() +
ggsave('plots/spatial_logcpue.png',g, width=7, height=4)
## png(paste0('plots/simulated_data', n_knots, '.png'), width=7, height=4, units='in', res=500)
## boxplot(df$logcpue.simulated, df$logcpue, names=c('Simulated', 'Observed'), ylab='CPUE')
## dev.off()

## Fit simple LM without space to visualize spatial residual patterns
mod.lm2 <- lm(logcpue.simulated~year+geartype+month+hooksize+depth, data=df)
## Run full spatio-temporal model
Version <- "models/spatiotemporal_cpue"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

## Make inputs
Data <- list(logcpue_i=df$logcpue.simulated,
             n_t=length(unique(df$year)),
             s_i=df.simulated$s_i-1,
             year_i=as.numeric(df$year)-1,
             month_i=as.numeric(df$month)-1,
             hooksize_i=as.numeric(df$hooksize)-1,
             geartype_i=as.numeric(df$geartype)-1,
             depth_i=df$depth,
             M0=spde$param.inla$M0, M1=spde$param.inla$M1,
             M2=spde$param.inla$M2)
Params <- list(intercept=4.4, beta_year=rep(0, length(levels(df$year))),
               beta_geartype=rep(0, length(levels(df$geartype))),
               beta_month=rep(0, length(levels(df$month))),
               beta_hooksize=rep(0, length(levels(df$hooksize))),
               beta_depth=0, ln_tau_O=log(1), ln_obs=log(1),
               ln_kappa=1, omega_s=rep(0,mesh$n),
               ln_tau_E=log(1),
               epsilon_st=matrix(0,nrow=nrow(Data$M0),ncol=Data$n_t))
Random <- c("omega_s", "epsilon_st")
# Build and run
Obj <- MakeADFun( data=Data, parameters=Params, random=Random )
Obj$env$beSilent()
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
             control=list(trace=1, eval.max=1e4, iter.max=1e4))
Opt$par
Opt[["final_diagnostics"]] <-
    data.frame( "Name"=names(Obj$par), "final_gradient"=Obj$gr(Opt$par))
Report <- Obj$report()
unlist( Report[c('Range','SigmaO','SigmaE')] )
## SD <- sdreport( Obj, bias.correct=TRUE)

df$resids.st <- Report$mu_i-df$logcpue.simulated
df$resids.lm <- fitted(mod.lm2)-df$logcpue.simulated
mean(df$resids.st)
mean(df$resids.lm)
range(c(df$resids.st, df$resids.lm))

df.temp <- subset(df, year==2001)
g <- ggplot(df.temp)+
    scale_color_continuous(low='blue', high='green', limits=c(-3,7)) +
     theme_bw() + facet_wrap('year')
g2 <- g + geom_point(aes(x=longitude, y=latitude, color=resids.lm), alpha=.5, size=.5)
ggsave('plots/map_resids_lm.png', g2, width=7, height=4)
g2 <- g + geom_point(aes(x=longitude, y=latitude, color=resids.st), alpha=.5, size=.5)
ggsave('plots/map_resids_st.png', g2, width=7, height=4)
