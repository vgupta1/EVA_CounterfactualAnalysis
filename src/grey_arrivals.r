######
# Generates predicted (counterfactual) time series of arrivals for
# greylisted countries had they NOT been greylisted
# Uses GBM and past arrivals data for predictions
# Includes estimated std errors for OPE estimates of value of greylisting
######

library(gbm)
library(ggplot2)
library(reshape2)

##########################
# read in files + preprocess
##########################

# read in SAMPLE full data
dat <- read.csv("../sample_data.csv")

# convert types, remove everything on/after 11/1 (end of CF analysis)
dat$date_entry <- as.Date(dat$date_entry, format = "%Y-%m-%d")
dat <- dat[dat$date_entry <= as.Date("2020-11-01", format = "%Y-%m-%d"),]
dat$to_test <- as.integer(as.character(dat$to_test) %in% c("true", "flag"))
dat$sent_for_test <- as.integer(!is.na(dat$sent_for_test) | dat$test_result %in% c("positive", "undefined"))

# subset vars of interest: country, arrivals & info for estimating show rate
arrivals <- dat[, c("country", "date_entry", "to_test", "sent_for_test")]
arrivals$count <- 1
arrivals <- aggregate(arrivals[,c("count", "to_test", "sent_for_test")],
                      by = list(arrivals$country, arrivals$date_entry), FUN = sum)
names(arrivals) <- c("country", "date", "plf", "to_test", "sent_for_test")

# identify countries as white/black/grey
ctry_white_list <- read.csv("https://raw.githubusercontent.com/ovelixasterix/greek_project/master/countries_allowed.csv", header = FALSE)
arrivals$white <- as.logical(as.character(arrivals$country) %in% ctry_white_list[,1])
grey_list_se <- read.csv("https://raw.githubusercontent.com/ovelixasterix/greek_project/master/grey_list_start_end.csv")
grey_list_se$start_date <- as.Date(grey_list_se$start_date, format = "%Y-%m-%d")
grey_list_se$end_date <- as.Date(grey_list_se$end_date, format = "%Y-%m-%d")
grey_list_se[is.na(grey_list_se$end_date),]$end_date <- as.Date("2022-01-01")
arrivals$grey <- FALSE
for(i in 1:nrow(arrivals)){
  ind <- which(grey_list_se$country == as.character(arrivals$country)[i])
  if(length(ind) > 0){
    arrivals$grey[i] <- as.logical(arrivals$date[i] >= grey_list_se[ind,]$start_date &
                                     arrivals$date[i] <= grey_list_se[ind,]$end_date)
  }
}

# smooth arrivals time series in a 7-day window, but only if same grey status
# smooth showup time series similarly
arrivals$sm_plf <- NA
arrivals$showup <- NA
for(i in 1:nrow(arrivals)){
  ind <- which(arrivals$country == arrivals$country[i] & arrivals$grey == arrivals$grey[i] &
                 abs(arrivals$date - arrivals$date[i]) <= 3)
  arrivals$sm_plf[i] <- mean(arrivals[ind,]$plf)
  arrivals$showup[i] <- mean(arrivals[ind,]$sent_for_test/arrivals[ind,]$to_test, na.rm = TRUE)
}

# add features on # arrivals in each category
arrivals$tot_white <- NA
arrivals$tot_black <- NA
for(i in 1:nrow(arrivals)){
  arrivals$tot_white[i] <- sum(arrivals[arrivals$date == arrivals$date[i] & arrivals$white,]$sm_plf)
  arrivals$tot_black[i] <- sum(arrivals[arrivals$date == arrivals$date[i] & !arrivals$white,]$sm_plf)
}
arrivals$white <- as.integer(arrivals$white)

arrivals$country <- as.factor(arrivals$country)

##########################
# predict greylist counterfactual arrivals
##########################

# test data
test <- arrivals[arrivals$grey,]
test$date <- as.integer(test$date - min(arrivals$date))

# train + validation set on non-grey
tmp <- arrivals[!arrivals$grey,]
tmp$date <- as.integer(tmp$date - min(arrivals$date))

# randomly split into train and validation set
n <- round(nrow(tmp)*.7)
ind <- sample(1:nrow(tmp))[1:n]
train <- tmp[ind,]
val <- tmp[-ind,]

# fit GBM model on train and predict on validation set
fit <- gbm(sm_plf ~ country + date + white + tot_white + tot_black, data = train, distribution = "gaussian", n.trees = 500, 
           bag.fraction = .75, cv.folds = 5, interaction.depth = 9)
pred_val <- predict(fit, newdata = val)
summary(fit)

# get standard deviation of residuals on validation set
residuals <- pred_val-val$sm_plf
hist(residuals)
sd_val <- sd(residuals)

# retrain on entire dataset (train + val)
fit <- gbm(sm_plf ~ country + date + white + tot_white + tot_black, data = rbind(train, val), distribution = "gaussian", n.trees = 500, 
           bag.fraction = .75, cv.folds = 5, interaction.depth = 9)

# test set predictions
pred_grey <- predict(fit, newdata = test)
test <- cbind(test, "pred_arr" = pred_grey, "sd_arr" = sd_val)

# predict on past data
grey <- unique(as.character(test$country))
dat <- arrivals[(as.character(arrivals$country) %in% grey) & !arrivals$grey,]
dat$date <- as.integer(dat$date - min(arrivals$date))
dat$pred_arr <- predict(fit, newdata = dat)
dat$sd_arr <- sd_val

# combine with test set
dat <- dat[,c("country", "date", "sm_plf", "pred_arr", "sd_arr", "showup")]
dat <- rbind(dat, test[,c("country", "date", "sm_plf", "pred_arr", "sd_arr", "showup")])
dat$date <- as.Date(dat$date + min(arrivals$date))

# scale arrivals, predictions & sd by showup rates
dat$showup <- pmin(dat$showup, 1)
dat$pred_arr <- dat$pred_arr * dat$showup
dat$sm_plf <- dat$sm_plf * dat$showup
dat$sd_arr <- dat$sd_arr * dat$showup

# make sure arrivals are non-negative
dat$pred_arr <- pmax(dat$pred_arr, 0)

# make sure counterfactual arrivals is below observed prevalence
# i.e., requiring negative-PCR pretest cannot increase arrivals
dat$pred_arr <- pmax(dat$pred_arr, dat$sm_plf)

write.csv(dat, "grey_predarr.csv")



