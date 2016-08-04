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
