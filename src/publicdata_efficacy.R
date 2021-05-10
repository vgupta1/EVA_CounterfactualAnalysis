######
# Try to predict private time series of EB estimates from public data
# Uses GBM and public data for predictions
######

library(gbm)
library(ggplot2)
library(ROCR)

##########################
# read in files + preprocess
##########################

# read in pre-processed public data file
X <- read.csv("../OtherData/Xall.csv", header = TRUE)

# remove row numbers (1st column)
X <- X[,-1]

# rename features & convert types
names(X) <- c("country", "date", paste("cases_minus", 1:20, sep=""), "cases_zero", 
              paste("cases_plus", 1:19, sep=""), paste("deaths_minus", 1:20, sep=""),
              "deaths_zero", paste("deaths_plus", 1:19, sep=""), paste("tests_minus", 1:20, sep=""),
              "tests_zero", paste("tests_plus", 1:19, sep=""))
X$date <- as.Date(X$date, format = "%Y-%m-%d")

# change NaNs to NAs
for(i in 3:ncol(X)){
  ind <- which(as.character(X[,i]) == "NaN")
  X[ind,i] <- NA
}

# identify countries as white/black, drop blacklisted
ctry_white_list <- read.csv("https://raw.githubusercontent.com/ovelixasterix/greek_project/master/countries_allowed.csv", header = FALSE)
X <- X[as.character(X$country) %in% ctry_white_list[,1],]

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
  # remove greylisted period
  ind <- which(tmp$date >= grey_list_se[i,]$start_date & tmp$date <= grey_list_se[i,]$end_date
               & as.character(tmp$country) == grey_list_se[i,]$country)
  if(length(ind) > 0) tmp <- tmp[-ind,]
  
  # remove * version of EB type for non-greylisted period
  ind <- which((tmp$date < grey_list_se[i,]$start_date | tmp$date > grey_list_se[i,]$end_date)
               & as.character(tmp$country) == grey_list_se[i,]$country
               & as.character(tmp$eb_type) == paste(grey_list_se[i,]$country, "*", sep=""))
  if(length(ind) > 0) tmp <- tmp[-ind,]
}

# merge with public data file
X <- merge(X, tmp, by=c("country", "date"))

# drop future columns
ind <- which(names(X) %in% c(paste("cases_plus", 1:19, sep=""), paste("deaths_plus", 1:19, sep=""), paste("tests_plus", 1:19, sep="")))
X <- X[,-ind]

# add features on sum of cases / deaths in last 14 days
ind <- which(names(X) %in% c("cases_zero", paste("cases_minus", 1:13, sep="")))
X$sum_cases_last14 <- rowSums(X[,ind])
ind <- which(names(X) %in% c("deaths_zero", paste("deaths_minus", 1:13, sep="")))
X$sum_deaths_last14 <- rowSums(X[,ind])
ind <- which(names(X) %in% c("tests_zero", paste("tests_minus", 1:13, sep="")))
X$sum_tests_last14 <- rowSums(X[,ind])
X$pos_rate_last14 <- X$sum_cases_last14/X$sum_tests_last14

X$country <- as.factor(X$country)

# split into train and test based on date, and convert date to integer after
n <- round(length(unique(X$date))*.8) + 1
trainEnd <- sort(unique(X$date))[n]
ind <- which(X$date < trainEnd)
trainX <- X[ind,]
testX <- X[-ind,]

##############################################
# fit GBM models on train and predict on test
##############################################

# Model 1: only sum of cases + deaths over last 14 days
feat <- c("eb_prev", "sum_cases_last14", "sum_deaths_last14")
fit1 <- gbm(eb_prev ~ ., data = trainX[,feat],
            distribution = "gaussian", n.trees = 500, 
            bag.fraction = .75, cv.folds = 5, interaction.depth = 9)
pred1 <- predict(fit1, newdata = testX[,feat])

# Model 2: sum of cases + deaths + tests + pos rate over last 14 days
feat <- c("eb_prev", "sum_cases_last14", "sum_deaths_last14", "sum_tests_last14", "pos_rate_last14")
fit2 <- gbm(eb_prev ~ ., data = trainX[,feat],
            distribution = "gaussian", n.trees = 500, 
           bag.fraction = .75, cv.folds = 5, interaction.depth = 9)
pred2 <- predict(fit2, newdata = testX[,feat])

# Model 3: cases + deaths time series
feat <- c("eb_prev", "cases_zero", paste("cases_minus", 1:13, sep=""), "sum_cases_last14",
          "deaths_zero", paste("deaths_minus", 1:13, sep=""), "sum_deaths_last14")
fit3 <- gbm(eb_prev ~ ., data = trainX[,feat],
            distribution = "gaussian", n.trees = 500, 
            bag.fraction = .75, cv.folds = 5, interaction.depth = 9)
pred3 <- predict(fit3, newdata = testX[,feat])

# Model 4: cases + deaths + tests time series
feat <- c("eb_prev", "cases_zero", paste("cases_minus", 1:13, sep=""), "sum_cases_last14",
          "deaths_zero", paste("deaths_minus", 1:13, sep=""), "sum_deaths_last14",
          "tests_zero", paste("tests_minus", 1:13, sep=""), "sum_tests_last14", "pos_rate_last14")
fit4 <- gbm(eb_prev ~ ., data = trainX[,feat],
            distribution = "gaussian", n.trees = 500, 
            bag.fraction = .75, cv.folds = 5, interaction.depth = 9)
pred4 <- predict(fit4, newdata = testX[,feat])

# Model 5: cases + deaths + tests time series & country-specific FE
feat <- c("eb_prev", "cases_zero", paste("cases_minus", 1:13, sep=""), "sum_cases_last14",
          "deaths_zero", paste("deaths_minus", 1:13, sep=""), "sum_deaths_last14",
          "tests_zero", paste("tests_minus", 1:13, sep=""), "sum_tests_last14", "pos_rate_last14",
          "country")
fit5 <- gbm(eb_prev ~ ., data = trainX[,feat],
            distribution = "gaussian", n.trees = 500, 
            bag.fraction = .75, cv.folds = 5, interaction.depth = 9)
pred5 <- predict(fit5, newdata = testX[,feat])

# get ROC curve around some cutoff
cutoff = 0.005 
testy_cut <- as.integer(testX$eb_prev >= cutoff)

colors = c("blue", "red", "forestgreen", "purple", "black")
auroc = c()
roc_results = data.frame("Model" = NA, "tpr" = NA, "fpr" = NA)
for(i in 1:5){
  tmp <- prediction(eval(parse(text=paste("pred", i, sep=""))), testy_cut)
  perf <- performance(tmp, measure = "auc")
  auroc <- c(auroc, perf@y.values[[1]])
  
  perf <- performance(tmp,"tpr","fpr")
  roc_results <- rbind(roc_results, data.frame("Model" = i,
                                               "fpr"= as.numeric(perf@x.values[[1]]),
                                               "tpr" = as.numeric(perf@y.values[[1]])))
  if(i==1){
    plot(perf, col="blue")
  } else {
    plot(perf, add=TRUE, col=colors[i])
  }
}

# viz / save ROC results
roc_results <- roc_results[-1,]
points(1:100/100, 1:100/100, type = "l", col = "black")
legend("bottomright", legend=c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"),
       col=colors, lty=c(1,1), cex=0.8)
# write.csv(roc_results, "../OPE_outputs/public_roc_results.csv", row.names = FALSE)

# viz / save AUROC results
tmp <- data.frame("Model" = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"),
                  "AUROC" = auroc)
ggplot(data=tmp, aes(x=Model, y=AUROC)) +
  geom_bar(stat="identity", fill="steelblue") + theme_minimal() +
  coord_cartesian(ylim=c(0,1))
# write.csv(auroc, "../OPE_outputs/public_auroc.csv", row.names = FALSE)

