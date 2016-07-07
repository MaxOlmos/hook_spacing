library(TMB)
library(plyr)

## Prep the data for JAGS:
## Vector of int representing which experiment it was from
hs <- read.csv("hs_data.csv")
hs <- subset(hs, used=="Y")             #
hs$hooks <- with(hs, hooks.per.skate*skates)
hs <- within(hs, {catch.per.hook=lbs/hooks; number.per.hook=number/hooks})
sum(hs$number)                   # seems to match the paper
sum(hs$lbs)                      # seems to match the paper
## A unique identifier for each row, for merging later
hs$rowID <- 1:nrow(hs)
## Normalize days
hs$date <- lubridate::mdy(as.character(hs$date))
temp <- ddply(hs, .(group), summarize,
              daynumber=as.numeric(difftime(date,min(date), units="days")),
              rowID=rowID)
hs <- merge(hs, temp, by="rowID", all.y=FALSE )
hs$group <- hs$group.x
temp <- hs$group.x <-  hs$group.y <- NULL
hs <- hs[with(hs, order(group, date)),]
hs$group <- as.numeric(hs$group)
## names(hs)
spacing <- hs$spacing
hs$group <- as.numeric(hs$group)
group <- hs$group
## The dependent variable (catch.per.hook for now) on log scale since it
## needs to be positive
## logy <- log(hs$catch.per.hook+.1)
yobs <- hs$catch.per.hook
log.yobs <- log(yobs+.1)
Nobs <- length(yobs)
Ngroup <- length(unique(group))
day <- hs$daynumber

g <- ggplot()+
    geom_line(data=hs, aes(x=spacing, y=catch.per.hook, group=daynumber,
              colour=daynumber), size=.5)  + xlim(0, 45)+ xlab("spacing (ft)")+
    geom_point(data=hs, aes(x=spacing, y=catch.per.hook, group=daynumber,
                  colour=daynumber), size=.75) + mytheme + ylab("catch/hook (lbs)")
g2 <- g + facet_wrap("ID", scales="fixed")
ggsave("../Plots/JAGS_Models/hs_data_fixed.png", g2, width=width, height=height)
g2 <- g + facet_wrap("ID", scales="free")
ggsave("../Plots/JAGS_Models/hs_data_free.png", g2, width=width, height=height)
## Same as above but on log scale
g <- ggplot()+
    geom_line(data=hs, aes(x=spacing, y=log(catch.per.hook+.1), group=daynumber,
              colour=daynumber), size=.5)  + xlim(0, 45)+ xlab("spacing (ft)")+
    geom_point(data=hs, aes(x=spacing, y=log(catch.per.hook+.1), group=daynumber,
                  colour=daynumber), size=.75) + mytheme + ylab("log(catch/hook)")
g2 <- g + facet_wrap("ID", scales="fixed")
ggsave("../Plots/JAGS_Models/hs_logdata_fixed.png", g2, width=width, height=height)
## g2 <- g + facet_wrap("ID", scales="free")
## ggsave("../Plots/JAGS_Models/hs_logdata_free.png", g2, width=width, height=height)
## Effect of daynumber
g <- ggplot(data=hs, aes(x=daynumber, y=log(.1+catch.per.hook)))+
    facet_wrap("ID", scales="fixed") +
    geom_point(size=1) + mytheme +
    ylab("log(catch/hook) (lbs)") +
    geom_line(aes(group=spacing), size=.5)  +
    stat_smooth(method="lm", fill="blue", colour="red",  size=.5, alpha=.5)+
    xlab("daynumber")
ggsave("../Plots/JAGS_Models/hs_data_dayeffect.png", g, width=width, height=height)


data.hs <-
    list(log_yobs=log.yobs, group=group-1, Nobs=Nobs, Ngroup=Ngroup, day=day,
         spacing=spacing)
Params <-
    list(logcpue_mean=0, logcpue_sd=1, sigma_obs_mean=-.5, sigma_obs_sd=.5,
         gamma=0, beta=.1,logcpue=rep(0, len=Ngroup), logsigma_obs=rep(0,
                                                        len=Ngroup))
Version <- "empirical_spacing"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
## Test model is working
Inputs <- list(Data=data.hs, Params=Params,
               random=NULL,#c('logsigma_obs', 'logcpue'),
               map=NULL)
lower <- rep(-Inf, len=length(unlist(Params)))
upper <- rep(Inf,  len=length(unlist(Params)))
names(upper) <- names(lower) <- names(unlist(Params))
lower['gamma'] <- 0; upper['gamma'] <- 1
lower['logcpue_sd'] <- 0; upper['logcpue_sd'] <- Inf
lower['sigma_obs_sd'] <- 0; upper['sigma_obs_sd'] <- Inf

Obj <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params,
                 random=Inputs$Random, map=Inputs$Map)
Obj$fn()
Obj$env$beSilent()
temp <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
               control=list(trace=10, eval.max=1e4, iter.max=1e4),
                            lower=lower, upper=upper)
test <- sdreport(Obj)
sd.df <- data.frame(spacing=1:70, value=test$value, sd=test$sd)
ggplot(sd.df, aes(spacing, value, ymax=value+2*sd, ymin=value-2*sd)) +
  geom_errorbar()
