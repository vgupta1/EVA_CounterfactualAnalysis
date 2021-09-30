############################################
# Script quantifies how much exploration Eva did under metric #1:
# how many tests were allocated to arms w/ eb_prev below current top number
# Note: test_alloc_per_country files in github are synthetic data for 2 dates only,
# results in paper are based on actual time series for entire summer
############################################

library(zoo)
library(ggplot2)

files <- list.files("OtherData/PastFits/")
files <- files[grep("test_alloc_per_country_", files)]

res <- data.frame("Date" = NA, "Explore" = NA, "Exploit" = NA)
res$Date <- as.Date(res$Date)
res <- res[-1,]

# loop through each file
for(i in 1:length(files)){
  # get date
  tmp <- strsplit(files[i], "_")[[1]][5]
  res[i,"Date"] <- as.Date(substr(tmp, 1, nchar(tmp) - 4))
  
  # read in file
  tmp <- paste("JonsFinalDataDump_26Nov/PastFits/", files[i], sep="")
  f <- read.csv(tmp)
  
  # get # of tests allocated
  nTests <- sum(f$numMarked)
  
  # get greedy allocations to top eb_prev types
  f <- f[order(f$eb_prev, decreasing = TRUE),]
  f$greedyMarked <- 0
  row = 1
  while(nTests > 0 & row <= nrow(f)){
    f$greedyMarked[row] = min(nTests, f$numArrived[row])
    nTests <- nTests - f$greedyMarked[row]
    row <- row + 1
  }
  
  # get # of tests allocated below min
  ind <- which(f$greedyMarked == 0)
  res[i, "Explore"] <- sum(f[ind,]$numMarked)
  res[i, "Exploit"] <- sum(f$numMarked) - sum(f[ind,]$numMarked)
  
}

# remove dates outside of consideration
res <- res[res$Date >= "2020-08-06" & res$Date <= "2020-11-01",]
sum(res$Explore)/(sum(res$Explore) + sum(res$Exploit))

# rolling avg, melt, restrict to peak season
res$Explore <- rollmean(res$Explore/(res$Explore + res$Exploit), k=3, fill=NA)
ggplot(res[c("Date", "Explore")], aes(x=Date,y=Explore)) + geom_line() +
  theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) + 
  ggtitle("Fraction of Tests Allocated to Exploration") + 
  scale_x_date(date_minor_breaks = "2 day")


# write.csv(res[c("Date", "Explore")], "explore_v1.csv", row.names=FALSE)
