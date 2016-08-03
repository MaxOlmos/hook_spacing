### The effect of hook spacing on longline catches.

### ------------------------------------------------------------
## Step 1: prepare workspace and load data
source('startup.R')
## source('data/load_data.R')

### ------------------------------------------------------------
## Step 2: Run the spatiotemporal model.
source('run_logbook.R')

### ------------------------------------------------------------
## Step 3: Run the Hamley & Skud model
source('run_empirical.R')

### ------------------------------------------------------------
## Step 4: Create plots, tables, and figures
source('make_figures.R')

### End of analysis
### ------------------------------------------------------------

## development code
x <- seq(.001, 70, len=1000)
ff <- function(x, beta, lambda) (1-exp(-beta*x)^lambda)/(1-exp(-beta*18)^lambda)
g <- function(x,beta,lambda){
  lines(x, ff(x, beta, lambda), type='l')
}
plot(18,1, ylim=c(0,2), xlim=range(x))
g(x,.05,10.1)




