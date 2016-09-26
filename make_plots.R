## This file loads the data, makes plots, tables and figures.

## Simulated trajectories with ignored hook spacing effect
fits.sim <- readRDS('results/fits.sim.RDS')
data.sim <- readRDS('results/data.sim.RDS')
densities.true <- ldply(1:3, function(i){
  x <- data.sim[[i]]
  data.frame(trend=c('original', 'flat', 'trend')[i], year=1996:2015, density=x$density_t)
  })
densities.sim <- ldply(fits.sim, function(x) x$sd.density)
x <- as.factor(densities.sim$.id)
levels(x) <- c('flat', 'original', 'trend')
densities.sim$trend <- x
densities.sim$density <- densities.sim$value
densities.sim$model <- "Predicted"
densities.true$model <- "Truth"
densities.both <- rbind(densities.true, densities.sim[, names(densities.true)])
densities.both <- ddply(densities.both, .(trend, model), mutate,
                        relative.density=density/mean(density),
                        test=mean(density))
g <- ggplot(densities.both, aes(year, relative.density, color=model,
                           group=model)) +
  geom_line() + facet_wrap('trend')
ggsave('plots/sim_trends.png', g, width=9, height=4)

## Make quick plots
fits.all <- readRDS('results/fits.all.RDS')
g <- plot.parameter.comparison(fits=fits.all[c(4,5,6)],
 level.name='model', levels=c('Nonparametric', 'Hamley & Skud', 'None'))
ggsave('plots/par_comparison_form.png', g, width=10, height=6)
g <- plot.cpue.comparison(fits.all)
ggsave('plots/cpue_comparison_form.png', g, width=10, height=6)
g <- plot.power.comparison(fits.all)
ggsave('plots/spacing_comparison.png', g, width=5, height=6)
g <- plot.resids.comparison(fits.all)
ggsave('plots/resids_comparison.png', g, width=5, height=6)


g <- ggplot(cpue.df, aes(year, cpue, ymin=lwr, ymax=upr, group=data, fill=data)) +
  geom_ribbon(alpha=.5)+facet_wrap('regarea', scales='free_y')
ggsave('plots/cpue_trends_vs_survey.png', g, width=ggwidth, height=ggheight)

g <- ggplot(df.summarized, aes(year, pct.catch, group=geartype, fill=geartype)) +
  geom_area() + ylab("Proportion Catch") + coord_flip()
ggsave('plots/proportion_geartype.png', g, width=ggwidth, height=ggheight)

g <- ggplot(spacing.sd, aes(spacing, value, ymax=upr, ymin=lwr)) + geom_ribbon() +
  facet_grid(model.name~form.name) + geom_line(col='red')
ggsave('plots/spacing_estimates.png', g, width=ggwidth, height=ggheight)
g <- ggplot(cpue.sd.normalized, aes(year, value2, ymax=upr2, ymin=lwr2)) + geom_ribbon() +
  facet_grid(model.name~form.name) + geom_line(col='red') + ylab('Normalized Abundance')
g <- g+ geom_hline(yintercept=0)
ggsave('plots/cpue_estimates.png', g, width=ggwidth, height=ggheight)

x1 <- droplevels(subset(spacing.sd, model=='ST' & form==2,
                 select=c(spacing, value, sd)))
x2 <- droplevels(subset(hs.results$uncertainty.df, select=c(spacing, value, sd)))
xx <- rbind.fill(cbind(model='Spatio-temporal',x1), cbind(model='Experimental', x2))
xx <- subset(xx, spacing < 43)
g <- ggplot(xx, aes(spacing, value, ymax=value+2*sd, ymin=value-2*sd,
                 fill=model, color=model)) + geom_ribbon(alpha=.5)
ggsave('plots/hs_estimates.png', g, width=ggwidth, height=ggheight)

xxxx


## The effect of resolution on the ST model
## knots <- c(50, 100, 150, 200, 300, 400, 500, 600, 700, 800)
knots <- c(50, 150, 200, 300, 400, 500, 600, 700, 800)
temp.list <- list(); k <- 1
for(n_knots in knots){
    x <- paste0('results/logbook.hs.', n_knots,'.RDS')
    if(file.exists(x)){
        y <- readRDS(x)
        temp <- with(y, data.frame(n_knots=n_knots, model=model, form=form, nll=report$nll_likelihood,
                    intercept=report$intercept, runtime=runtime))
        temp.list[[k]] <- cbind(temp,
        rbind(data.frame(par='runtime', value=y$runtime, sd=0), y$sd.par[1:5,]));
        k <- k+1
    }
}
res <- do.call(rbind, temp.list)
g <- ggplot(res, aes(n_knots, value, ymin=value-2*sd, ymax=value+2*sd)) +
    geom_line() + geom_errorbar() + facet_wrap('par', scales='free')
ggsave('plots/resolution_effect.png', g, width=7, height=5)

## Load main results
empirical.results <- readRDS('results/empirical.results.RDS')
logbook.re.results <- readRDS('results/logbook.re.results.RDS')
logbook.hs.results <- readRDS('results/logbook.hs.results.RDS')

## Make CSV files for tables
sd.pars <- rbind(cbind(model='empirical', empirical.results$sd.df),
                 cbind(model='logbook.re', logbook.re.results$sd.par),
                 cbind(model='logbook.hs', logbook.hs.results$sd.par))
sd.pars <- subset(sd.pars, par !='cph_t')
sd.pars$table <- with(sd.pars,
 paste0(round(value,3), ' (', round(value-1.96*sd,3), '-', round(value+1.96*sd,3), ')'))
g <- ggplot(sd.pars, aes(model, value, ymin=value-1.96*sd, ymax=value+1.96*sd)) + geom_pointrange() + facet_wrap('par')
ggsave('plots/parameter_estimates_by_model.png', g, width=7, height=5)
write.table(x=sd.pars, file='results/sd.table.csv', sep=',', row.names=FALSE)

## Plot effect hook comparisons for the three models
uncertainty.df <- rbind.fill(data.frame(model='empirical', empirical.results$uncertainty.df),
      data.frame(model='logbook.re', logbook.re.results$sd.spacing),
      data.frame(model='logbook.hs', logbook.hs.results$sd.spacing))
hs.original <- data.frame(spacing=1:max(uncertainty.df$spacing))
hs.original$value=(1-exp(-.06*hs.original$spacing))/(1-exp(-.06*18))
hs.original$sd <- 0
g <- ggplot(uncertainty.df, aes(spacing, value, ymax=value+2*sd, ymin=value-2*sd)) +
  geom_errorbar() + geom_line(lwd=2) + facet_wrap('model') +
      geom_line(data=hs.original, aes(x=spacing, y=value), col='red') +
      ylab('Effective Hook') + xlab("Hook Spacing (ft)") #+ xlim(0,50)
ggsave('plots/effective_hook_comparisons.png', g, width=7, height=3)

## QQ plots of residuals for the logbook models to see if they're fitting
## well
xx <- data.frame(model='logbook.re', qqnorm(logbook.re.results$report$resids, plot=FALSE))
yy <- data.frame(model='logbook.hs', qqnorm(logbook.hs.results$report$resids, plot=FALSE))
g <- ggplot(rbind(xx,yy), aes(x,y)) +geom_point()+ facet_wrap('model') +
    geom_abline(slope=1, intercept=0)
ggsave('plots/qq_plots.png', g, width=7, height=3)


## ## Make some exploratory plots of the logbook data
## g <- ggplot(df, aes(year, cph)) + geom_violin()
## ggsave('plots/raw_cph.png', g, width=ggwidth, height=ggheight)
## g <- ggplot(df, aes(year, log(cph))) + geom_violin()
## ggsave('plots/raw_logcpue.png', g, width=ggwidth, height=ggheight)
## g <- ggplot(df, aes(spacing)) + geom_bar() + facet_wrap('geartype')
## ggsave('plots/spacings.png', g, width=ggwidth, height=ggheight)
## df.geartype <- ddply(df, .(year, geartype), summarize,
##                      count=length(geartype))
## df.geartype <- ddply(df.geartype, .(year), mutate, pct=count/sum(count))
## g <- ggplot(df.geartype, aes(year,pct, fill=geartype)) + geom_bar(stat='identity')
## ggsave('plots/geartype_trends.png', g, width=ggwidth, height=ggheight)

## ## Exploratory plots of the emprical data
## g <- ggplot()+
##     geom_line(data=hs, aes(x=spacing, y=catch.per.hook, group=daynumber,
##               colour=daynumber), size=.5)  + xlim(0, 45)+ xlab("spacing (ft)")+
##     geom_point(data=hs, aes(x=spacing, y=catch.per.hook, group=daynumber,
##                   colour=daynumber), size=.75) + mytheme + ylab("catch/hook (lbs)")
## g2 <- g + facet_wrap("ID", scales="fixed")
## ggsave("../Plots/JAGS_Models/hs_data_fixed.png", g2, width=width, height=height)
## g2 <- g + facet_wrap("ID", scales="free")
## ggsave("../Plots/JAGS_Models/hs_data_free.png", g2, width=width, height=height)
## ## Same as above but on log scale
## g <- ggplot()+
##     geom_line(data=hs, aes(x=spacing, y=log(catch.per.hook+.1), group=daynumber,
##               colour=daynumber), size=.5)  + xlim(0, 45)+ xlab("spacing (ft)")+
##     geom_point(data=hs, aes(x=spacing, y=log(catch.per.hook+.1), group=daynumber,
##                   colour=daynumber), size=.75) + mytheme + ylab("log(catch/hook)")
## g2 <- g + facet_wrap("ID", scales="fixed")
## ggsave("../Plots/JAGS_Models/hs_logdata_fixed.png", g2, width=width, height=height)
## ## g2 <- g + facet_wrap("ID", scales="free")
## ## ggsave("../Plots/JAGS_Models/hs_logdata_free.png", g2, width=width, height=height)
## ## Effect of daynumber
## g <- ggplot(data=hs, aes(x=daynumber, y=log(.1+catch.per.hook)))+
##     facet_wrap("ID", scales="fixed") +
##     geom_point(size=1) + mytheme +
##     ylab("log(catch/hook) (lbs)") +
##     geom_line(aes(group=spacing), size=.5)  +
##     stat_smooth(method="lm", fill="blue", colour="red",  size=.5, alpha=.5)+
##     xlab("daynumber")
## ggsave("../Plots/JAGS_Models/hs_data_dayeffect.png", g, width=width, height=height)


## ## Process the results and make plots
## Report.list <- readRDS('results/report_models.RDS')
## Results <- do.call(rbind, lapply(Report.list, function(x)
##   data.frame(x[c('intercept', 'beta_depth','SigmaE', 'Range', 'Sigma',
##                  'SigmaO', 'jnll', 'time')])))
## Results$model <- c('M1', "M2", "M3")
## SD.list <- readRDS('sd_models.RDS')
## cpue.trends <- (lapply(SD.list, function(x)
##   data.frame(year=1995:2012, cpue=x['value'], se=x['sd'])))
## cpue.trends <- melt(cpue.trends, id.vars=c('year', 'value', 'sd'))
## names(cpue.trends) <- c('year', 'cpue', 'sd', 'model')
## cpue.trends <- within(cpue.trends, {
##                         lwr=cpue-sd
##                         upr=cpue+sd})

## ## CPUE trends with and without errors
## g <- ggplot(cpue.trends, aes(year, cpue, group=model, color=model)) +
##   geom_line() + ylab("Relative CPUE")
## ggsave('plots/cpue_trends.png', g, width=ggwidth, height=ggheight)
## g <- ggplot(cpue.trends, aes(year, cpue, ymin=lwr, ymax=upr)) +
##   geom_linerange(color=gray(.5)) + geom_point()+ylab("Relative CPUE") + facet_wrap('model')
## ggsave('plots/cpue_trends_errorbars.png', g, width=ggwidth, height=ggheight)

## ## Model parameters and meta data
## Results.long <- melt(Results, id.vars='model')
## g <- ggplot(Results.long, aes(model, (value), group=variable)) +
##   facet_wrap('variable', scales='free_y') + geom_line()
## ggsave('plots/results_by_model.png', g, width=ggwidth, height=ggheight)
## qqplots <- ldply(Report.list, function(x){
##                    xx <- qqnorm(y=x$resids, plot=FALSE)
##                  return(data.frame(xx))})
## g <- ggplot(qqplots, aes(x,y)) + geom_point() + facet_wrap('.id') + geom_abline(slope=1,intercept=0)
## ggsave('plots/qqresids.png', g, width=ggwidth, height=ggheight)

## ## Plots of resids to see if there's a spatial trend removed
## df$resids.M1 <- Report.list[['M1']]$resids
## df$resids.M2 <- Report.list[['M2']]$resids
## df$resids.M3 <- Report.list[['M3']]$resids
## col.scale <-  scale_colour_distiller(palette='Spectral', limits=c(-3,3)) #scale_color_continuous(low='blue', high='red')
## myxlim <- xlim(-155,-145); myylim <- ylim(56,61);
## ## Some global ones
## g <- ggplot(df, aes(longitude, latitude, color=resids.M1)) +
##   geom_point(alpha=1, size=.15) + col.scale
## ggsave('plots/spatial_errors_M1.png',g, width=ggwidth, height=ggheight)
## g <- ggplot(df, aes(longitude, latitude, color=resids.M2)) +
##   geom_point(alpha=1, size=.15) + col.scale
## ggsave('plots/spatial_errors_M2.png',g, width=ggwidth, height=ggheight)
## g <- ggplot(df, aes(longitude, latitude, color=resids.M3)) +
##   geom_point(alpha=1, size=.15) + col.scale
## ggsave('plots/spatial_errors_M3.png',g, width=ggwidth, height=ggheight)

## for(yr in unique(df$year)){
## g <- ggplot(subset(df, year==yr), aes(longitude, latitude, color=resids.M1)) +
##   geom_point(alpha=1, size=.5) + col.scale + myxlim+myylim
## ggsave(paste0('plots/annual/spatial_errors_',yr,'_M1.png'),g, width=7, height=4)
## g <- ggplot(subset(df, year==yr), aes(longitude, latitude, color=resids.M2)) +
##   geom_point(alpha=1, size=.5) + col.scale + myxlim+myylim
## ggsave(paste0('plots/annual/spatial_errors_',yr,'_M2.png'),g, width=7, height=4)
## g <- ggplot(subset(df, year==yr), aes(longitude, latitude, color=resids.M3)) +
##   geom_point(alpha=1, size=.5) + col.scale + myxlim+myylim
## ggsave(paste0('plots/annual/spatial_errors_',yr,'_M3.png'),g, width=7, height=4)
## }
