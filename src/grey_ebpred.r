######
# Generates predicted (counterfactual) time series of EB estimates for
# greylisted countries had they NOT been greylisted
# Uses GBM and public data for predictions
# Includes estimated std errors for OPE estimates of value of greylisting
######

library(gbm)
library(ggplot2)
library(reshape2)


##########################
# read in files + preprocess
##########################

# read in pre-processed public data file
X <- read.csv("../OtherData/Xall.csv", header = TRUE)

# remove row numbers (1st column)
X <- X[,-1]

# rename features & convert date to correct type
names(X) <- c("country", "date", paste("cases_minus", 1:20, sep=""), "cases_zero", 
              paste("cases_plus", 1:19, sep=""), paste("deaths_minus", 1:20, sep=""),
              "deaths_zero", paste("deaths_plus", 1:19, sep=""), paste("tests_minus", 1:20, sep=""),
              "tests_zero", paste("tests_plus", 1:19, sep=""))
X$date <- as.Date(X$date, format = "%Y-%m-%d")

# identify countries as white/black
# note that all greylisted countries were white
ctry_white_list <- read.csv("https://raw.githubusercontent.com/ovelixasterix/greek_project/master/countries_allowed.csv", header = FALSE)
X$white <- as.integer(as.character(X$country) %in% ctry_white_list[,1])

# identify when countries were greylisted
grey_list_se <- read.csv("https://raw.githubusercontent.com/ovelixasterix/greek_project/master/grey_list_start_end.csv")
grey_list_se$start_date <- as.Date(grey_list_se$start_date, format = "%Y-%m-%d")
grey_list_se$end_date <- as.Date(grey_list_se$end_date, format = "%Y-%m-%d")
grey_list_se[is.na(grey_list_se$end_date),]$end_date <- as.Date("2022-01-01")
X$grey <- FALSE
for(i in 1:nrow(X)){
  ind <- which(grey_list_se$country == as.character(X$country)[i])
  if(length(ind) > 0){
    X$grey[i] <- as.logical(X$date[i] >= grey_list_se[ind,]$start_date &
                              X$date[i] <= grey_list_se[ind,]$end_date)
  }
}

# merge with OPE outputs file of private EB predictions
tmp <- read.csv("../OPE_outputs/hist_eb_timeseries_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv")
tmp <- tmp[c("country", "eb_type", "date", "eb_prev")]
tmp$date <- as.Date(tmp$date, format = "%Y-%m-%d")
for(i in 1:nrow(grey_list_se)){
  # remove non-* version of EB type for greylisted period
  ind <- which(tmp$date >= grey_list_se[i,]$start_date & tmp$date <= grey_list_se[i,]$end_date
               & as.character(tmp$country) == grey_list_se[i,]$country
               & as.character(tmp$eb_type) == grey_list_se[i,]$country)
  if(length(ind) > 0) tmp <- tmp[-ind,]
  
  # remove * version of EB type for non-greylisted period
  ind <- which((tmp$date < grey_list_se[i,]$start_date | tmp$date > grey_list_se[i,]$end_date)
               & as.character(tmp$country) == grey_list_se[i,]$country
               & as.character(tmp$eb_type) == paste(grey_list_se[i,]$country, "*", sep=""))
  if(length(ind) > 0) tmp <- tmp[-ind,]
}

# merge with public data file
X2 <- merge(X, tmp, by=c("country", "date"))

X2$country <- as.factor(X2$country)
X2$eb_type <- NULL

##########################
# predict greylist counterfactual EB prevalence
##########################

# test set, convert date to a numeric since some start date
X_grey <- X2[X2$grey,]
X_grey$grey <- NULL
X_grey$date <- as.integer(X_grey$date - min(X2$date))

# train + validation set on non-grey
tmp <- X2[!X2$grey,]
tmp$grey <- NULL
tmp$date <- as.integer(tmp$date - min(X2$date))

# randomly split into train and validation set
n <- round(nrow(tmp)*.7)
ind <- sample(1:nrow(tmp))[1:n]
trainX <- tmp[ind, -which(names(tmp) %in% c("eb_prev"))]
trainy <- tmp[ind, "eb_prev"]
valX <- tmp[-ind, -which(names(tmp) %in% c("eb_prev"))]
valy <- tmp[-ind, "eb_prev"]

# fit GBM model on train and predict on validation set
fit <- gbm(trainy ~ ., data = trainX, distribution = "gaussian", n.trees = 500, 
           bag.fraction = .75, cv.folds = 5, interaction.depth = 9)
pred_val <- predict(fit, newdata = valX)
summary(fit)

# get standard deviation of residuals on validation set
residuals <- pred_val-valy
hist(residuals)
sd_val <- sd(residuals)

# retrain on entire dataset (train + val)
fit <- gbm(c(trainy, valy) ~ ., data = rbind(trainX, valX), distribution = "gaussian", n.trees = 500, 
           bag.fraction = .75, cv.folds = 5, interaction.depth = 9)

# test set predictions
pred_grey <- predict(fit, newdata = X_grey[,-which(names(X_grey) %in% c("eb_prev"))])
X_grey <- cbind(X_grey, "pred_prev" = pred_grey, "sd_prev" = sd_val)

# predict on past data
grey <- unique(as.character(X_grey$country))
dat <- X2[(as.character(X2$country) %in% grey) & !X2$grey,]
dat$date <- as.integer(dat$date - min(X2$date))
dat$pred_prev <- predict(fit, newdata = dat[,-which(names(dat) %in% c("grey", "eb_prev"))])
dat$sd_prev <- sd_val

# combine with test set
dat <- dat[,c("country", "date", "eb_prev", "pred_prev", "sd_prev")]
dat <- rbind(dat, X_grey[,c("country", "date", "eb_prev", "pred_prev", "sd_prev")])
dat$date <- as.Date(dat$date + min(X2$date))

# make sure predicted prevalence is non-negative
dat$pred_prev <- pmax(dat$pred_prev, 0)

# make sure counterfactual prevalence is above observed prevalence
# i.e., requiring negative-PCR pretest cannot increase prevalence
dat$pred_prev <- pmax(dat$pred_prev, dat$eb_prev)

write.csv(dat, "../OPE_outputs/grey_eb_preds.csv")

