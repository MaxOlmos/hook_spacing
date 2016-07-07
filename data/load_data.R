
## ldlog <- read.csv("ldlog.csv")
## save(ldlog, file="ldlog.Rdata")
print("This file loads ldlog and cleans it up")
load("ldlog.Rdata")
effective.skates <- function(hooks.per.skate, hook.spacing, skates=1){
    (skates*hooks.per.skate)*1.52*(1-exp(-.06*hook.spacing))/100
}

## Make some global changes to the data.
## Clean up the area codes
ldlog$regcde <- as.character(ldlog$regcde)
ldlog$regcde <- gsub(pattern=" ", replacement="", x=ldlog$regcde)
ldlog$regcde <- gsub(pattern="4EE|CLS", replacement="4E", x=ldlog$regcde)
ldlog$regcde <- gsub(pattern="4C|4D|4E", replacement="4CDE", x=ldlog$regcde)
ldlog$regcde[is.na(ldlog$regcde)] <- "Unknown"
ldlog$regcde <- factor(ldlog$regcde, levels=sort(unique(ldlog$regcde)))
ldlog$statarea <- as.factor(as.character(ldlog$statarea))
## Clean up gear codes. Make a new column for snap, fixed , or other gear
## type and collapse some of the obscure ones into an "other" category.
ldlog$grcde <- as.character(ldlog$grcde)
ldlog$grcde[is.na(ldlog$grcde)] <- "Unknown"
ldlog$grcde <- gsub(pattern="BC|CH|FS|PT|TP|TR|UC|UL", replacement="Other", x=ldlog$grcde)
ldlog$grcde <- factor(ldlog$grcde, levels=sort(unique(ldlog$grcde)))
ldlog$geartype <- as.character(ldlog$grcde)
ldlog$geartype <- gsub("SN|SS|SU|SG|SK", "Snap", ldlog$geartype)
ldlog$geartype <- gsub("FH", "Fixed", ldlog$geartype)
ldlog$geartype <- gsub("AU", "Autoline", ldlog$geartype)
ldlog$geartype <- with(ldlog, factor(geartype, sort(unique(geartype))))
## Hook size comes in two different scales. Using a quick merge to make a
## new column where the size is consistent. From Tom: ".  Note there are
## two scales, one that goes ¡§low to high¡¨ and another that goes ¡§high to
## low¡¨.  Size 3=16, 4=15, 5=14, 6=13¡K ". Assuming for now that the system
## switches at 8.
ldlog$hksiz.original <- ldlog$hksiz
for(i in 1:8) {
    ldlog$hksiz[ldlog$hksiz.original==(1:8)[i]] <- (18:11)[i]
}
## ----------
## Effective skates.
## # hooks can be calculated from 2 of the 3 variables. Currently they
## filter out those that have an implied hook spacing more than 5%
## different than what is reported. I repeated that here using the
## "valid.effort" column which is TRUE/FALSE
## mean(is.na(ldlog$hkskt))
## spacing.na <- is.na(ldlog$hkspc) | ldlog$hkspc == 0
## length.na <- is.na(ldlog$sktlen) | ldlog$sktlen == 0
## hooksperskate.na <- is.na(ldlog$hkskt) | ldlog$hkskt == 0
## table(spacing.na,  hooksperskate.na)
## Check that the implied hook spacing works out to be the same or close
## enough
hkskt_implied <- with(ldlog, sktlen/hkspc)
ldlog$valid.effort <- with(ldlog, abs(hkskt_implied-hkskt) < .05*pmin(hkskt_implied, hkskt))
## with(ldlog, plot(hkskt, abs(hkskt_implied-hkskt),  col=valid.effort+1, pch=".",))
## abline(h=0)
## Now calculate the actual number of hooks and effective skates
ldlog$effskate <- with(ldlog,
   ifelse(valid.effort, effective.skates(hkskt, hkspc, skates=skthld), NA))
ldlog$cpue.original <- with(ldlog, catwgt/effskthld)
ldlog$cpue <- with(ldlog, catwgt/effskate)
ldlog$cpue[is.infinite(ldlog$cpue)] <- NA
ldlog$cpue.corrected <- with(ldlog, ifelse(geartype=="Snap", cpue*1.35, cpue))
## Some crazy outliers in CPUE which make the plots realy hard to read, so
## I created a truncated version where high values were capped at 500.
ldlog$truncated.cpue <- with(ldlog, pmin(cpue, 500))
## catch per hook (cph)
ldlog$cph <- with(ldlog, catwgt/(skthld*hkskt)) # catch/hook=weight/(skates*hooks/skate)
## ldlog$lost.gear <- with(ldlog, ifelse(is.na(sktlst), 0, sktlst/(sktlst+skthld)))
## Clean up the dates
ldlog$catdate <- ymd(as.character(ldlog$catdate))
ldlog$month <- as.factor(month(ldlog$catdate))
## Clean up the locations, from Tom:
## 'Lat/Lat2/Lon/Lon2 ¡V coordinates in degrees decimal-minutes
## (eg. 542186=54¢X21.86¡¦) Note: prior to 2011, coords were rounded to the
## nearest minute.'
ldlog$lon1num <- -with(ldlog, as.numeric(substr(lon, 1, 3)) +
                       as.numeric(substr(lon, 4, 7))/100/60)
ldlog$lon2num <- -with(ldlog, as.numeric(substr(lon2, 1, 3)) +
                       as.numeric(substr(lon2, 4, 7))/100/60)
ldlog$lat1num <- with(ldlog, as.numeric(substr(lat, 1, 2)) +
                       as.numeric(substr(lat, 3, 6))/100/60)
ldlog$lat2num <- with(ldlog, as.numeric(substr(lat2, 1, 2)) +
                       as.numeric(substr(lat2, 3, 6))/100/60)
## Longitude values beyond 180 are stored in a strange way. From Ian: "I
## talked with Aaron, our survey DB guy, and he said the longitudes past
## 180 degrees just go back down with 100 added to them.  To put them back
## on a linear scale I think you should be able to do something like:
## Plot_long = (180 - (DB_long - 100))+180". His formula is not quite
## right, but I think this does it:
ldlog$lon1num <- with(ldlog, ifelse(lon1num < -180, -360-(lon1num+100),
                                    lon1num))
ldlog$lon2num <- with(ldlog, ifelse(lon2num < -180, -360-(lon2num+100),
                                    lon2num))
## Combine together into a single one for plotting. Note that the duplicate
## columns give the starting and ending locations of the hual. To combine
## together I'll just take the simple average of these numbers. If one one
## is given (the other NA) I'll take that one.
ldlog$longitude <-  rowMeans(ldlog[, c("lon1num", "lon2num")], na.=TRUE)
ldlog$longitude[is.nan(ldlog$longitude)] <- NA
ldlog$latitude <-  rowMeans(ldlog[, c("lat1num", "lat2num")], na.=TRUE)
ldlog$latitude[is.nan(ldlog$latitude)] <- NA
## plot(latitude~longitude, data=ldlog)
ldlog$haul.length <-
    with(ldlog, sqrt( (lon1num-lon2num)^2 + (lat1num-lat2num)^2))
## Clean up the reason codes
ldlog$target <- ifelse(is.na(ldlog$rsncde), "Halibut", ldlog$rsncde)
ldlog$target <- gsub("3|5|6|7|8|9", "Mixed", x=ldlog$target)
ldlog$target <- gsub("2|4", "Abnormal/Incomplete", x=ldlog$target)
ldlog$year <- ldlog$logyr
ldlog.annual <- ddply(ldlog, .(year,regcde), summarize,
                     pct.coords=100*(1-mean(is.na(longitude) | is.na(latitude))),
                     total.sets=length(year))
ldlog.annual.wide <- melt(ldlog.annual, id.vars=c('year', 'regcde'))
g <- ggplot(ldlog.annual.wide, aes(year, value, color=regcde, group=regcde)) + geom_line() +
    facet_grid(variable~., scales='free')
ggsave('plots/pct.coords.png', g, width=ggwidth, height=ggheight)

### end of processing whole data set

### For project just using certain area, and need to futher cleanup the data.
df <- subset(ldlog, year >= 1995 & year < 2013 & regcde=='3A' & longitude
             >-156 & longitude < -135 & latitude > 56 & target=='Halibut' &
             geartype != "Unknown" & cpue.corrected>0)
df$depth <- rowMeans(df[,c('dep1', 'dep2')], na.rm=TRUE)
df$hooksize <- factor(df$hksiz)
df$year <- factor(df$year)
df$cpue <- df$cpue.corrected            # using version corrected for snap gear
df$logcpue <- log(df$cpue)
df$spacing <- df$hkspc
df <- subset(df, select=c(cph, cpue, logcpue, year, spacing, month, geartype, statarea,
                   hooksize, depth, longitude, latitude))
df <- data.frame(droplevels(na.omit(df)))
g <- ggplot(df, aes(geartype, spacing)) + geom_violin()
ggsave('plots/spacings.png', g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(longitude, latitude, color=depth)) + geom_point(alpha=.1, size=.01)
ggsave('plots/spatial_depth.png', g, width=ggwidth, height=ggheight)
g <- ggplot(df, aes(longitude, latitude, color=statarea)) + geom_point(alpha=.1, size=.1)
ggsave('plots/spatial_statarea.png', g, width=ggwidth, height=ggheight)

message("Saving data to data.RDS for later use")
saveRDS(df, file="data.RDS")
## rm(df, ldlog, ldlog.annual, g, i, effective.skates, hkskt_implied)
## End of processing data
## ------------------------------------------------------------
