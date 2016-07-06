### The effect of hook spacing.

### ------------------------------------------------------------
## Step 1: prepare workspace and explore data
source('startup.R')
## source('load_data.R')
df <- readRDS(file='data.RDS')
str(df)
n_years <- length(unique(df$year))
df$spacing <- round(df$spacing)
## Make some exploratory plots of the data
g <- ggplot(df, aes(year, cph)) + geom_violin()
ggsave('plots/raw_cph.png', g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(year, logcpue)) + geom_violin()
ggsave('plots/raw_logcpue.png', g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(spacing)) + geom_bar() + facet_wrap('geartype')
ggsave('plots/spacings.png', g, width=ggwidth, height=ggheight)

### ------------------------------------------------------------
## Step 2: Run the spatiotemporal model.
source('run_cpue.R')

### ------------------------------------------------------------
## Step 3: Run the Hamley & Skud model
source('run_hs.R')

### ------------------------------------------------------------
## Step 4: Create plots, tables, and figures



### End of analysis
### ------------------------------------------------------------
