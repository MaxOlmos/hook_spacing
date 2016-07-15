## load and process the empirical data
hs <- read.csv("data/hs_data.csv")
hs <- subset(hs, used=="Y")             #
hs$hooks <- with(hs, hooks.per.skate*skates)
hs <- within(hs, {catch.per.hook=lbs/hooks; number.per.hook=number/hooks})
## A unique identifier for each row, for merging later
hs$rowID <- 1:nrow(hs)
## Normalize days
hs$date <- lubridate::mdy(as.character(hs$date))
hs <- ddply(hs, .(group), mutate,
              daynumber=as.numeric(difftime(date,min(date), units="days")))
hs <- hs[with(hs, order(group, date)),]
hs$group <- as.numeric(hs$group)
## The dependent variable (catch.per.hook for now) on log scale since it
## needs to be positive. For now adding .1 but should revisit later

data.hs <-
    list(log_yobs=log(hs$catch.per.hook+.1), group=hs$group-1,
         Nobs=nrow(hs), Ngroup=length(unique(hs$group)), day=hs$daynumber,
         spacing=hs$spacing)
Params <-
    list(logcpue_mean=0, logcpue_sd=1, sigma_obs_mean=-.5, sigma_obs_sd=.5,
         gamma=0, beta=.1, lambda=1, logcpue=rep(0, len=data.hs$Ngroup), logsigma_obs=rep(0,
                                                        len=data.hs$Ngroup))
Map <- list(lambda=factor(NA))
Version <- "models/empirical_spacing"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
## Set bounds for parameters via limits in optimizer
lower <- rep(-Inf, len=length(unlist(Params)))
upper <- rep(Inf,  len=length(unlist(Params)))
names(upper) <- names(lower) <- names(unlist(Params))
lower['gamma'] <- 0; upper['gamma'] <- 1
lower['lambda'] <- 0; upper['lambda'] <- 5
lower['logcpue_sd'] <- 0; upper['logcpue_sd'] <- Inf
lower['sigma_obs_sd'] <- 0; upper['sigma_obs_sd'] <- Inf
Obj <- MakeADFun(data=data.hs, parameters=Params,
                 random=c('logsigma_obs', 'logcpue'), map=Map)
Obj$env$beSilent()
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
               control=list(trace=10, eval.max=1e4, iter.max=1e4),
                            lower=lower, upper=upper)
temp <- sdreport(Obj)
xx <- data.frame(par=names(temp$value), value=temp$value, sd=temp$sd)
uncertainty.df <- subset(xx, par=='ehook')
uncertainty.df$spacing <- 1:nrow(uncertainty.df)
sd.df <- subset(xx, par!= 'ehook')
empirical.results <- list(Obj=Obj, Opt=Opt, sdreport=temp,
                          uncertainty.df=uncertainty.df, sd.df=sd.df)
saveRDS(empirical.results, file='results/empirical.results.RDS')

## Clean up
dyn.unload(dynlib(Version))
rm(Obj, Opt, temp, uncertainty.df, empirical.results, lower, upper,
   Version, Params)
