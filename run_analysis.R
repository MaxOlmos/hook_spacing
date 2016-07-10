### The effect of hook spacing on longline catches.

### ------------------------------------------------------------
## Step 1: prepare workspace and load data
source('startup.R')
## source('load_data.R')
df <- readRDS(file='data/data.RDS')
str(df)
n_years <- length(unique(df$year))
df$spacing <- round(df$spacing)

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
