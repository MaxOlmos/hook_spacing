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
Inputs <- make.inputs.experimental(hs)
Version <- "models/empirical_spacing"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

## Make object and optimize
Obj <- with(Inputs, MakeADFun(data=Data, parameters=Params,
                 random=Random, map=Map))
Obj$par
Obj$env$beSilent()
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
               control=list(trace=10, eval.max=1e4, iter.max=1e4),
                            lower=Inputs$lower, upper=Inputs$upper)
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
rm(Obj, Opt, temp, uncertainty.df, empirical.results, Inputs)
