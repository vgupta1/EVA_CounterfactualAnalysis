rm(list=ls())
setwd("~/Dropbox (Penn)/Greek Targeted Testing (Kimon-Vishal)/CounterFactualAnalysis/")

dat <- read.csv("grey_eb_preds.csv")
dat$X <- NULL
dat$date <- as.Date(dat$date, format = "%Y-%m-%d")

dat2 <- read.csv("grey_arrivs.csv")
dat2$X <- NULL
dat2$date <- as.Date(dat2$date, format = "%Y-%m-%d")

grey_list_se <- read.csv("https://raw.githubusercontent.com/ovelixasterix/greek_project/master/grey_list_start_end.csv")
grey_list_se$start_date <- as.Date(grey_list_se$start_date, format = "%Y-%m-%d")

lag = 9 # alter if desired

for(i in 1:nrow(grey_list_se)){
  # irrelevant dates
  ind <- which(dat$country == grey_list_se$country[i] & (dat$date <= grey_list_se$start_date[i]
               | dat$date > (grey_list_se$start_date[i] + lag)))
  dat <- dat[-ind,]
  ind <- which(dat2$country == grey_list_se$country[i] & (dat2$date <= grey_list_se$start_date[i]
                                                         | dat2$date > (grey_list_se$start_date[i] + lag)))
  dat2 <- dat2[-ind,]
}

# merge
df <- merge(dat[c("country", "date", "pred_prev", "sd_prev")],
            dat2[c("country", "date", "pred_arr", "sd_arr")],
            by = c("country", "date"))
df$peak = as.integer(df$date < as.Date("2020-10-01"))
ind <- which(df$peak == 1)

# numerically simulate sum, drawing arr and eb predictions from independent normals
B = 1000000
res = matrix(0, nrow = nrow(df), ncol = B)
for(i in 1:nrow(df)){
  res[i,] = rnorm(B, mean = df[i,"pred_prev"], sd = df[i,"sd_prev"]) * rnorm(B, mean = df[i,"pred_arr"], sd = df[i,"sd_arr"])
}

# peak
res2 <- colSums(res[ind,])
sd(res2)

# off-peak
res2 <- colSums(res[-ind,])
sd(res2)


