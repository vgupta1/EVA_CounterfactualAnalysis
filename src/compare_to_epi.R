############################################
# Script evaluates how many infections we would have caught had we relied
# on importance sampling countries by a public data metric (cases, deaths or positivity rate)
# and used that to allocate tests at each point of entry on each date
# uses IPW instead of direct estimate
# assigns global probs to each type based on public data
############################################

library(reshape2)
library(ggplot2)
library(zoo)

##############
## READ IN FILES, CLEAN PUBLIC DATA METRICS
##############

# read in OPE info
ope <- read.csv("OPE_outputs/ope_dat_6Dec_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv")
ope$estNumArrivals <- round(ope$estNumArrivals)
ope$date_entry <- as.Date(ope$date_entry)

# restrict to dates of interest
ope <- ope[ope$date_entry >= "2020-08-06" & ope$date_entry <= "2020-11-01",]

# read in preprocessed public data #s
pub <- read.csv("Other Data/Xall_new.csv")
pub$date <- as.Date(pub$date)

# restrict to dates of interest
pub <- pub[pub$date >= "2020-08-06" & pub$date <= "2020-11-01",]

# get yesterday's case rate (if unavailable get last available one)
tmp <- pub[,rev(grep("cases\\.", names(pub)))]
pub$cases <- as.vector(apply(tmp, 1, function(x) x[!is.na(x) & !is.infinite(x)][1]))

# get yesterday's positivity rate (if unavailable get last available one)
tmp <- pub[,rev(grep("cases\\.", names(pub)))]/pub[,rev(grep("tests\\.", names(pub)))]
pub$pos <- as.vector(apply(tmp, 1, function(x) x[!is.na(x) & !is.infinite(x)][1]))

tmp <- pub[,rev(grep("deaths\\.", names(pub)))]
pub$deaths <- as.vector(apply(tmp, 1, function(x) x[!is.na(x) & !is.infinite(x)][1]))

# remove other columns
pub <- pub[,c("country", "date", "cases", "deaths", "pos")]

# for each metric, zero our negatives & impute if value missing/infinite
for(i in 3:5){
  ind <- which(pub[,i] < 0)
  pub[ind,i] <- 0
  
  ind <- which(is.na(pub[,i]) | is.infinite(pub[,i]))
  if(length(ind) > 0){
    pub[ind,i] <- mean(pub[-ind,i])
  }
}

publicPerf <- function(metrics, ope, pub){
  pub <- pub[c("country","date", metrics)]
  names(pub)[3] <- "metric"
  
  ##############
  ## COMPUTE PROPENSITIES FOR PUBLIC DATA
  ##############
  
  dates <- unique(ope$date_entry)
  ctries <- unique(ope$country)
  
  # store results: propensity score (tau_1 in notes) for each (country, date) pair
  res <- data.frame(matrix(0, ncol = length(ctries) + 1, nrow = length(dates)))
  names(res) <- c("date", ctries)
  res$date <- dates
  
  arrival_matrix <- res
  xk_matrix <- res
    
  # loop through dates
  for(i in 1:length(dates)){
    
    # get date
    dt <- dates[i]
    
    # get pass manifest for that date (country, entry, numarrivals)
    tmp <- ope[ope$date_entry == dt, c("country", "point_entry", "estNumArrivals")]
    pass <- aggregate(.~country+point_entry, tmp, sum)
    pass <- pass[pass$estNumArrivals > 0,] # remove types w/ 0 arrivals
    
    # get test budget for every point of entry
    tmp <- ope[ope$date_entry == dt, c("point_entry", "numTestsPerformed")]
    budget <- aggregate(.~point_entry, tmp, sum)
    budget <- budget[budget$numTestsPerformed > 0,] # remove ports w/ 0 tests
    
    # get public case-based risk for that day
    pub2 <- pub[pub$date == dt, c("country", "metric")]
    
    # impute values for countries w/ arrivals but no public data
    ind <- which(!pass$country %in% pub2$country)
    if(length(ind) > 0){
      tmp <- data.frame("country" = unique(pass$country[ind]))
      tmp$metric = mean(pub2$metric)
      pub2 <- rbind(pub2, tmp)
    }

    # loop through countries
    for(j in 1:length(ctries)){
      
      pass2 <- merge(pass, pub2)
      budget2 <- merge(budget, pass2[pass2$country == ctries[j],],
                       by="point_entry")
      
      if(nrow(budget2) == 0){ # no arrivals from this country
        res[i,j+1] <- NA
      }
      else{
        for(e in 1:nrow(budget2)){
          ake <- budget2$estNumArrivals[e]
          ak <- sum(budget2$estNumArrivals)
          Te <- budget2$numTestsPerformed[e]
          
          tmp <- pass2[pass2$point_entry == budget2$point_entry[e],]
          normx <- sum(tmp$estNumArrivals * tmp$metric)
          xk <- tmp[tmp$country == ctries[j],]$metric
          
          val <- min(Te*xk/normx, 1) # ensure probs cant exceed 1
          res[i,j+1] <- res[i,j+1] + ake*val/ak 
          
          # print(ake*Te*xk/(ak * normx))
        }
        arrival_matrix[i,j+1] <- ak
        xk_matrix[i,j+1] <- xk
      }
    }
  }

  # melt & write out prop scores
  tmp <- melt(res,id="date")
  f <- paste("prob_public_", metrics, "_test", sep="")
  names(tmp) <- c("date", "country", f)
  ind <- which(is.na(tmp[,3]))
  tmp <- tmp[-ind,]
  write.csv(tmp, paste("Other Data/", f, ".csv", sep=""), row.names=FALSE)
  
  # melt and write out xk values
  xk_matrix <- melt(xk_matrix,id="date")
  f <- paste(metrics, "_metric", sep="")
  names(xk_matrix) <- c("date", "country", f)
  write.csv(xk_matrix, paste("Other Data/", f, ".csv", sep=""), row.names=FALSE)
  
  # get frac of tests performed relative to eva (scaling factor)
  ind <- which(res$date <= "2020-10-01")
  nTests <- sum(res[ind,-1]*arrival_matrix[ind,-1], na.rm = TRUE)
  scalingPeak <- sum(ope[ope$date_entry <= "2020-10-01",]$numTestsPerformed)/nTests

  ind <- which(res$date > "2020-10-01")
  nTests <- sum(res[ind,-1]*arrival_matrix[ind,-1], na.rm = TRUE)
  scalingOffPeak <- sum(ope[ope$date_entry > "2020-10-01",]$numTestsPerformed)/nTests
  
  ##############
  ## COMPUTE # CASES CAUGHT FOR PUBLIC DATA USING IPW
  ##############
  
  # get gittins testing probs
  ope2 <- aggregate(.~ date_entry + country + prob_gittins_test,
                    ope[c("date_entry", "country", "prob_gittins_test", "numPositivesFound")], sum)
  names(ope2)[1] <- "date"
  
  # merge with public policy testing probs & compute propensity scores
  ope2 <- merge(ope2, tmp, by=c("date", "country"))
  ope2$prop_score <- ope2[,5]/ope2$prob_gittins_test
  ind <- which(is.na(ope2$prop_score) | is.infinite(ope2$prop_score))
  ope2 <- ope2[-ind,]
  
  # remove scores exceeding the 97.5% quantile
  ind <- which(ope2$prop_score >= quantile(ope2$prop_score, probs=.975))
  ope2 <- ope2[-ind,]
  
  # compute num caught and aggregate over time
  ope2$numCaught <- ope2$prop_score * ope2$numPositivesFound
  numCaught <- aggregate(.~date, ope2[c("date", "numCaught", "numPositivesFound")], sum)
  names(numCaught) <- c("date", "numPublicCaught", "numGittinsCaught")

  return(list(numCaught, c(scalingPeak, scalingOffPeak)))
}

######
# run fn for each metric
######

# CASES
tmp <- publicPerf(c("cases"), ope, pub)
numCaught <- tmp[[1]]
scalings <- tmp[[2]]

# scale by tests performed in peak/offpeak
ind <- which(numCaught$date <= "2020-10-01")
numCaught$numPublicCaught[ind] <- numCaught$numPublicCaught[ind]*scalings[1]
numCaught$numPublicCaught[-ind] <- numCaught$numPublicCaught[-ind]*scalings[2]

sum(numCaught[ind,]$numGittinsCaught)/sum(numCaught[ind,]$numPublicCaught)
sum(numCaught[-ind,]$numGittinsCaught)/sum(numCaught[-ind,]$numPublicCaught)

names(numCaught)[2:3] <- c("Cases", "Eva")

# DEATHS
tmp <- publicPerf(c("deaths"), ope, pub)
tmpCaught <- tmp[[1]]
scalings <- tmp[[2]]

# scale by tests performed in peak/offpeak
ind <- which(tmpCaught$date <= "2020-10-01")
tmpCaught$numPublicCaught[ind] <- tmpCaught$numPublicCaught[ind]*scalings[1]
tmpCaught$numPublicCaught[-ind] <- tmpCaught$numPublicCaught[-ind]*scalings[2]

sum(tmpCaught[ind,]$numGittinsCaught)/sum(tmpCaught[ind,]$numPublicCaught)
sum(tmpCaught[-ind,]$numGittinsCaught)/sum(tmpCaught[-ind,]$numPublicCaught)

names(tmpCaught)[2] <- "Deaths"
tmpCaught$numGittinsCaught <- NULL
numCaught <- merge(numCaught, tmpCaught)

# POSITIVITY
tmp <- publicPerf(c("pos"), ope, pub)
tmpCaught <- tmp[[1]]
scalings <- tmp[[2]]

# scale by tests performed in peak/offpeak
ind <- which(tmpCaught$date <= "2020-10-01")
tmpCaught$numPublicCaught[ind] <- tmpCaught$numPublicCaught[ind]*scalings[1]
tmpCaught$numPublicCaught[-ind] <- tmpCaught$numPublicCaught[-ind]*scalings[2]

sum(tmpCaught[ind,]$numGittinsCaught)/sum(tmpCaught[ind,]$numPublicCaught)
sum(tmpCaught[-ind,]$numGittinsCaught)/sum(tmpCaught[-ind,]$numPublicCaught)

names(tmpCaught)[2] <- "Positivity"
tmpCaught$numGittinsCaught <- NULL
numCaught <- merge(numCaught, tmpCaught)

# TOTAL SUMMER: improvement of EVA relative to what public data would have caught
sum(numCaught$Eva)/colSums(numCaught[,c("Cases", "Deaths", "Positivity")])

# PEAK
tmp <- numCaught[numCaught$date <= "2020-10-01",]
sum(tmp$Eva)/colSums(tmp[,c("Cases", "Deaths", "Positivity")])

# OFFPEAK
tmp <- numCaught[numCaught$date > "2020-10-01",]
sum(tmp$Eva)/colSums(tmp[,c("Cases", "Deaths", "Positivity")])

## VISUALIZE

# normalize by Eva's performance, take rolling avg
melt <- numCaught
melt$Cases <- rollmean(melt$Eva/melt$Cases, k=7, fill=NA)
melt$Deaths <- rollmean(melt$Eva/melt$Deaths, k=7, fill=NA)
melt$Positivity <- rollmean(melt$Eva/melt$Positivity, k=7, fill=NA)
melt$Eva <- NULL

# rolling avg, melt, restrict to peak season
melt <- melt(melt,id="date")
melt <- melt[melt$date <= "2020-10-01",]
ggplot(melt, aes(x=date,y=value,colour=variable,group=variable)) + geom_line() +
  ylab("Relative Infections Caught by Eva") + xlab("Date") + 
  scale_x_date(date_minor_breaks = "2 day")






