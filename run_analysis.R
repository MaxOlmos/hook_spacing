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
source('run_experimental.R')

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


## explore gamma model to make sure I understand it

k <- 5
theta <- 50/k
xmax <- qgamma(p=.999, shape=k, scale=theta)
x <- seq(0.001, xmax, len=1000)
par(mfrow=c(1,2))
plot(x, dgamma(x, shape=k, scale=theta))
abline(v=theta*k)
theta <- 5000/k
xmax <- qgamma(p=.999, shape=k, scale=theta)
x <- seq(qgamma(0.001, k, scale=theta), xmax, len=1000)
plot(x, dgamma(x, shape=k, scale=theta))
abline(v=theta*k)

N <- 5000
CV <- .01
beta0 <- 3
beta1 <- 2
x <- runif(N, 0, 500)
nu <- x*beta1+beta0
y <- rgamma(N, shape=1/CV^2, scale=nu*CV^2)
nll <- function(par){
  CV=par[1]
  pred=x*par[3]+par[2]
  -sum(dgamma(x=y, shape=1/CV^2, scale=pred*CV^2, log=TRUE))
}
opt <- optim(par=c(CV, beta0, beta1), fn=nll, lower=c(0,0,0), method='L-BFGS-B')
## calculate Pearson residuals
pred=x*opt$par[3]+opt$par[2]
resids <- (pred-y)
var <- (1/CV^2)*(pred*CV^2)^2
resids.pearson <- (pred-y)/(pred*CV)
par(mfrow=c(1,2))
qqnorm(resids)
qqnorm(resids.pearson)

## matches wikipedia?
x <- seq(0.001, 20, len=10000)
mat <- cbind(c(1,2,3,5,9,7.5,.5), c(2,2,2,1,.5,1,1))
yy <- ldply(1:nrow(mat), function(xx)
  cbind(x=x, row=xx, shape=mat[xx,1], theta=mat[xx,2], density=dgamma(x, shape=mat[xx,1], scale=mat[xx,2])))
ggplot(yy, aes(x, density, group=row, color=row)) + geom_line() + ylim(0,.5)

