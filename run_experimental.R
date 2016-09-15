source("startup.R")

## load and process the empirical data
hs <- read.csv("data/hs_data.csv")
hs <- droplevels(subset(hs, used=="Y" & lbs>0))
hs$hooks <- with(hs, hooks.per.skate*skates)
hs$catch <- hs$lbs
## Normalize days
hs$date <- lubridate::mdy(as.character(hs$date))
hs <- ddply(hs, .(group), mutate,
              daynumber=as.numeric(difftime(date,min(date), units="days")))
hs <- hs[with(hs, order(group, date)),]
hs$group <- as.numeric(hs$group)

## Prepare TMB inputs
Data <-
  list(n_i=nrow(hs), n_s=length(unique(hs$group)),
       catch_i=hs$catch, hooks_i=hs$hooks, day_i=hs$daynumber,
       spacing_i=hs$spacing, site_i=hs$group-1)
Params <-
    list(eta_mean=1, eta_sd=1, sigma_mean=-.5, sigma_sd=.5,
         gamma=0.1, beta=.5, lambda=1,
         eta_s=rep(1, len=Data$n_s), sigma_s=rep(1,len=Data$n_s))
Map <- list(lambda=factor(NA))          # can't estimate this so fix at 1
Version <- "models/empirical_spacing"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
## Set bounds for parameters via limits in optimizer
lower <- rep(-Inf, len=length(unlist(Params)))
upper <- rep(Inf,  len=length(unlist(Params)))
names(upper) <- names(lower) <- names(unlist(Params))
lower['gamma'] <- 0; upper['gamma'] <- 1
lower['lambda'] <- 0; upper['lambda'] <- 5
lower['eta_sd'] <- 0; upper['eta_sd'] <- Inf
lower['sigma_sd'] <- 0; upper['sigma_sd'] <- Inf

## Make object and optimize
Obj <- MakeADFun(data=Data, parameters=Params,
                 random=c('eta_s', 'sigma_s'), map=Map)
Obj$env$beSilent()
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
               control=list(trace=10, eval.max=1e4, iter.max=1e4),
                            lower=lower, upper=upper)
temp <- sdreport(Obj)
xx <- data.frame(par=names(temp$value), value=temp$value, sd=temp$sd)
uncertainty.df <- subset(xx, par=='hook_power')
uncertainty.df$spacing <- 1:nrow(uncertainty.df)
sd.df <- subset(xx, par!= 'hook_power')
ggplot(uncertainty.df, aes(spacing, value, ymax=value+2*sd, ymin=value-2*sd)) + geom_pointrange()
empirical.results <- list(Obj=Obj, Opt=Opt, sdreport=temp,
                          uncertainty.df=uncertainty.df, sd.df=sd.df)
saveRDS(empirical.results, file='results/empirical.results.RDS')

## Clean up
dyn.unload(dynlib(Version))
rm(Obj, Opt, temp, uncertainty.df, empirical.results, lower, upper,
   Version, Params)
