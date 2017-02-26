## This file loads and preps the data for making figures plots and tables

## logbook results for the stacked poly plot
df.summarized <- readRDS('results/df.summarized.RDS')

survey <- readRDS('data/survey.RDS')
fits.areas <- readRDS('results/fits.areas.RDS')
cpue.areas <- ldply(seq_along(regareas), function(x)
  cbind(data='logbook', regarea=regareas[x], fits.areas[[x]]$sd.density))
cpue.areas$se <- cpue.areas$sd
cpue.areas <- ddply(cpue.areas, .(regarea), mutate,
      cpue=value/mean(value),
      lwr=lwr/mean(value),
      upr=upr/mean(value))
cpue.areas <- cpue.areas[, names(survey)]
cpue.df <- rbind.fill(survey, cpue.areas)

fits.all <- readRDS('results/fits_form_vs_model.RDS')
par.sd.table <- ldply(fits.all, function(x)
  cbind(form=x$form, model=x$model, x$sd.par))
par.sd.table <- droplevels(subset(par.sd.table, model %in% c('NS', 'ST')))

options(scipen=6)
par.sd.table$table <-
  with(par.sd.table,
       paste0(formatC(value,  width=3, digits=3), " (", formatC(sd, width=3,digits=3),  ")"))
par.sd.table$table
## par.sd.table$table2 <-
##   with(par.sd.table,
##        paste0(formatC(value, digits=3, sci=TRUE), " (", formatC(lwr,
##                                                               digits=3, sci=TRUE), "-",
##               formatC(upr, digits=3, sci=TRUE), ")"))
## par.sd.table$table <- with(par.sd.table,
##   ifelse((value!=0 & abs(value)<.01) | is.na(value), table2, table1))

## quick table of parameter estimates and SDs
par.sd.table.wide <- dcast(par.sd.table, par2~model+form, value.var='table')
write.table(x=par.sd.table.wide, 'results/table_logbook.csv',
            row.names=FALSE, sep=",")

spacing.sd <- ldply(fits.all, function(x){
  cbind(form=x$form, form.name=form.name(x$form), model=x$model,
  model.name=model.name(x$model),  x$sd.spacing)
  })
spacing.sd <- droplevels(subset(spacing.sd, model %in% c('NS', 'ST') & form!=3))

cpue.sd <- ldply(fits.all, function(x){
  cbind(form=x$form, form.name=form.name(x$form), model=x$model,
  model.name=model.name(x$model), x$sd.density)})

cpue.sd <- droplevels(subset(cpue.sd, model %in% c('NS', 'ST') & form!=1))
cpue.sd.normalized <- ddply(cpue.sd, .(model, form), mutate,
                            value2=value/mean(value),
                            upr2=upr/mean(value),
                            lwr2=lwr/mean(value))

## experimental results
hs.results <- readRDS('results/experimental.results.RDS')
hs.table <- hs.results$sd.df
hs.table$table <- with(hs.table,
       paste0(formatC(value,  width=3, digits=3), " (", formatC(sd, width=3,digits=3),  ")"))
hs.table <- hs.table[, c('par', 'table')]
write.table(x=hs.table, 'results/table_experimental.csv',
            row.names=FALSE, sep=",")
