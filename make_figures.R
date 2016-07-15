## This file loads the data, makes plots, tables and figures.

## The effect of resolution on the ST model
knots <- c(50, 100, 150, 200, 300, 400, 500, 600, 700, 800)
knots <- c(50, 150, 200, 300, 400, 500, 600, 700, 800)
temp.list <- list(); k <- 1
for(n_knots in knots){
    x <- paste0('results/logbook.hs.', n_knots,'.RDS')
    if(file.exists(x)){
        y <- readRDS(x)
        temp <- with(y,
                     data.frame(n_knots=n_knots, model=model, form=form, nll=report$nll_likelihood,
                                intercept=report$intercept, runtime=runtime))
        temp.list[[k]] <- cbind(temp, y$sd.par[1:5,]); k <- k+1
    }
}
res <- do.call(rbind, temp.list)
ggplot(res, aes(n_knots, value, ymin=value-2*sd, ymax=value+2*sd)) +
    geom_line() + geom_errorbar() + facet_wrap('par', scales='free')
ggplot(res, aes(n_knots, runtime)) + facet_wrap('par', scales='free') + geom_line()

empirical.results <- readRDS('results/empirical.results.RDS')
logbook.re.results <- readRDS('results/logbook.re.results.RDS')
logbook.hs.results <- readRDS('results/logbook.hs.results.RDS')

## Plot effect hook comparisons for the three models
uncertainty.df <- rbind.fill(data.frame(model='empirical', empirical.results$uncertainty.df),
      data.frame(model='logbook.re', logbook.re.results$uncertainty.df),
      data.frame(model='logbook.hs', logbook.hs.results$uncertainty.df))
hs.original <- data.frame(spacing=1:max(uncertainty.df$spacing))
hs.original$value=(1-exp(-.06*hs.original$spacing))/(1-exp(-.06*18))
hs.original$sd <- 0
g <- ggplot(uncertainty.df, aes(spacing, value, ymax=value+2*sd, ymin=value-2*sd)) +
  geom_errorbar() + geom_line(lwd=2) + facet_wrap('model') +
      geom_line(data=hs.original, aes(x=spacing, y=value), col='red') +
      ylab('Effective Hook') + xlab("Hook Spacing (ft)") + xlim(0,50)
ggsave('plots/effective_hook_comparisons.png', g, width=7, height=5)

## QQ plots of residuals for the logbook models to see if they're fitting
## well
xx <- data.frame(model='logbook.re', qqnorm(logbook.re.results$report$resids, plot=FALSE))
yy <- data.frame(model='logbook.hs', qqnorm(logbook.hs.results$report$resids, plot=FALSE))
ggplot(rbind(xx,yy), aes(x,y)) +geom_point()+ facet_wrap('model') +
    geom_abline(slope=1, intercept=0)

## Make some exploratory plots of the logbook data
g <- ggplot(df, aes(year, cph)) + geom_violin()
ggsave('plots/raw_cph.png', g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(year, log(cph))) + geom_violin()
ggsave('plots/raw_logcpue.png', g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(spacing)) + geom_bar() + facet_wrap('geartype')
ggsave('plots/spacings.png', g, width=ggwidth, height=ggheight)
df.geartype <- ddply(df, .(year, geartype), summarize,
                     count=length(geartype))
df.geartype <- ddply(df.geartype, .(year), mutate, pct=count/sum(count))
g <- ggplot(df.geartype, aes(year,pct, fill=geartype)) + geom_bar(stat='identity')
ggsave('plots/geartype_trends.png', g, width=ggwidth, height=ggheight)

## Exploratory plots of the emprical data
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


## Process the results and make plots
Report.list <- readRDS('results/report_models.RDS')
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
