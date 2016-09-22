### this file contains intermediate tests I ran during development. It
### probably won't run straight through but the code chunks might be useful
### later so I'm saving them here.

## Explore difference in the likelihoods
NS11 <- run.logbook(n_knots=10, model='NS' , form=1, likelihood=1)
NS21 <- run.logbook(n_knots=10, model='NS' , form=2, likelihood=1)
NS12 <- run.logbook(n_knots=10, model='NS' , form=1, likelihood=2)
NS22 <- run.logbook(n_knots=10, model='NS' , form=2, likelihood=2)
plot.spacing.fit(list(NS11, NS21, NS12, NS22), TRUE)
par(mfrow=c(1,2))
qqnorm(NS21$report$resids)
qqnorm(NS22$report$resids)

S11 <- run.logbook(n_knots=400, model='S' , form=1, likelihood=1)
S21 <- run.logbook(n_knots=400, model='S' , form=2, likelihood=1)
S12 <- run.logbook(n_knots=400, model='S' , form=1, likelihood=2)
S22 <- run.logbook(n_knots=400, model='S' , form=2, likelihood=2)
plot.spacing.fit(list(S11, S21, S12, S22), TRUE)
par(mfrow=c(1,2))
qqnorm(S21$report$resids)
qqnorm(S22$report$resids)

S1 <- run.logbook(n_knots=200, model='S' , form=1, likelihood=1)
S2 <- run.logbook(n_knots=200, model='S' , form=2, likelihood=1)
ST1 <- run.logbook(n_knots=200, model='ST' , form=1, likelihood=1)
ST2 <- run.logbook(n_knots=200, model='ST' , form=2, likelihood=1)
plot.spacing.fit(list(NS1, NS2, S1, S2, ST1, ST2), TRUE)

## Increase spatial resolution to see when it's sufficiently large
for(n_knots in c(50, 100, 150, 200, 300, 400, 500, 600, 700, 800)){
    form <- 2
    print(n_knots)
    temp <- run.logbook(n_knots=n_knots, model='ST', form=form, trace=0)
    if(form==1) saveRDS(temp, file=paste0('results/logbook.re.', n_knots,'.RDS'))
    if(form==2) saveRDS(temp, file=paste0('results/logbook.hs.', n_knots,'.RDS'))
}


## Likelihoods
knots <- 1000
mlog <- run.logbook(n_knots=knots, model='NS', form=3, vessel=TRUE)
mgam <- run.logbook(n_knots=knots, model='NS', form=3, vessel=TRUE, likelihood=2)
minvg <- run.logbook(n_knots=knots, model='NS', form=3, vessel=TRUE, likelihood=3)
fits <- list(mlog, mgam, minvg)
temp <- ldply(1:length(fits), function(x)
  data.frame(likelihood=fits[[x]]$likelihood, form=fits[[x]]$form,
             resids=fits[[x]]$report$resids, preds=fits[[x]]$report$preds))
temp$likelihood <- factor(temp$likelihood, labels=c('Lognormal', 'Gamma', 'Inverse Gaussian'))
g <- ggplot(temp, aes(preds, resids)) + geom_point(alpha=.1, size=.1) +
  xlab('Predicted Catch') + ylab('Pearson Residual') + facet_wrap('likelihood') +
  scale_x_log10() + ylim(-5,5) + theme_bw() + geom_hline(yintercept=0)
ggsave('plots/outliers_normal.png', g, width=map.width, height=map.height)
df.temp <- df
df.temp$resids=mlog$report$resids
g <- ggplot(droplevels(subset(df.temp, abs(resids)>4)),
            aes(longitude, latitude, color=geartype)) +
  geom_point(alpha=.7) + ggtitle('Lognormal model: Pearson residuals > 4')
rm(df.temp)
ggsave('plots/map_outliers_normal.png', g, width=map.width, height=map.height)
g <- plot.parameter.comparison(fits,
     level.name='likelihood', levels=c('Lognormal', 'Gamma', 'Inv. Gaussian'))
ggsave('plots/par_comparison_likelihood.png', g, width=ggwidth, height=ggheight)

## Vessel effect
knots <- 1000
v0 <- run.logbook(n_knots=knots, model='ST', form=2, vessel_effect=FALSE, trace=10)
v1 <- run.logbook(n_knots=knots, model='ST', form=2, vessel_effect=TRUE, trace=10)
plot.parameter.comparison(list(v0, v1),
     level.name='vessel', levels=c('No Vessels', 'Vessel Effect'))
ggsave('plots/par_comparison_vessel.png')
png('plots/vessel_qqplots.png', res=500, units='in', width=14, height=8)
par(mfrow=c(1,2))
with(v0$report, {qqnorm(resids, main='No vessel effect'); qqline(resids)})
with(v1$report, {qqnorm(resids, main='Vessel effect'); qqline(resids)})
dev.off()

## Data filtering
knots <- 1000
## data.unfiltered <- readRDS(file='data/data_unfiltered.RDS')
## data.unfiltered$spacing <- round(data.unfiltered$spacing)
d1 <- run.logbook(n_knots=knots, model='ST', form=2, vessel=FALSE)
df.temp <- df
df <- df.unfiltered
d2 <- run.logbook(n_knots=knots, model='ST', form=2, vessel=FALSE)
g <- plot.parameter.comparison(list(d1,d2),
     level.name='data', levels=c('Filtered', 'Unfiltered'))
ggsave('plots/par_comparison_data.png', g, width=8, height=6)
df <- df.temp
rm(df.unfiltered, df.temp)

## Simulate some simple data without effects for anything except spacing
## to make sure the model recovers it
beta <- 0.08
lambda <- .1
max_ehook <- 1/(1-exp(-18*beta)^lambda)
df.sim <- df
df.sim$spacing <- sample(1:40, size=nrow(df.sim), replace=TRUE)
df.sim$catch <- with(df.sim, hooks*(1-exp(-beta*spacing)^lambda)*
      max_ehook*exp(.5)*exp(rnorm(nrow(df.sim), 0,.8)))
knots <- 10
df <- df.sim
sim1 <- run.logbook(n_knots=knots, model='NS', form=1, vessel=FALSE)
sim2 <- run.logbook(n_knots=knots, model='NS', form=2, vessel=FALSE)
g <- plot.parameter.comparison(list(sim1, sim2),
     level.name='form', levels=c('Smoother', 'H&S', 'None'))
g
g <- plot.power.comparison(fits=list(sim1, sim2))
g



## Spatial effect
knots <- 1000
ns <- run.logbook(n_knots=knots, model='NS', form=2, vessel=TRUE)
s <- run.logbook(n_knots=knots, model='S', form=2, vessel=TRUE)
st <- run.logbook(n_knots=knots, model='ST', form=2, vessel=TRUE)
fits.model <- list(ns,s,st)
saveRDS(fits.model, file='results/fits.model.RDS')
fits.model <- readRDS('results/fits.model.RDS')
g <- plot.parameter.comparison(fits.model,
     level.name='model', levels=c('NS', 'S', 'ST'))
ggsave('plots/par_comparison_model.png')

## Spacing effect for NS model
knots <- 10
ns1 <- run.logbook(n_knots=knots, model='NS', form=1, vessel=FALSE)
ns2 <- run.logbook(n_knots=knots, model='NS', form=2, vessel=FALSE)
ns3 <- run.logbook(n_knots=knots, model='NS', form=3, vessel=FALSE)
fits.spacing <- list(ns1,ns2,ns3)
saveRDS(fits.spacing, file='results/fits.spacing.RDS')
fits.spacing <- readRDS('results/fits.spacing.RDS')
g <- plot.parameter.comparison(fits.spacing,
     level.name='form', levels=c('Smoother', 'H&S', 'None'))
ggsave('plots/par_comparison_model.png', g, width=12, height=6)
g <- plot.spacing.comparison(fits.spacing)
ggsave('plots/spacing_comparison.png', g, width=5, height=6)
g <- plot.power.comparison(fits.spacing)
ggsave('plots/power_comparison.png', g, width=5, height=6)
## plot(df$spacing, ns2$report$resids, pch='.', col=rgb(0,0,0,.1))
