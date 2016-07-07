## This file loads the data, makes plots, tables and figures.

empirical.results <- readRDS('results/empirical.results.RDS')

ggplot(empirical.results$uncertainty.df, aes(spacing, value, ymax=value+2*sd, ymin=value-2*sd)) +
  geom_errorbar() + geom_line(lwd=2)
