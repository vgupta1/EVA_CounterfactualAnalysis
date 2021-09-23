#### 
#  New Variance Specifications for Revision
#  Implements new variance spec and a cleaner OPE analysis for public policy stuff. 
###

library(tidyverse)
library(lubridate)
library(modelr)

PROP_CUT_OFF = 34.7
END_OF_SEASON = ymd("2020-10-01")

# dat should be like ope_dat spec with certain extra columns
#   #prob_counterfactual_test = tau_1 for the public test.
#   #metric = x_k(t) for relevant type/time
# q_prop_thresh the quantile of propensity scores. (.975)  above this are dropped.

compute_var_auc <- function(dat, q_prop_thresh = .975, DELTA_AC = 8, test_scalings = NA){
  if (is.na(test_scalings)){
    test_scalings <- tribble(~is_peak, ~scaling, 
                           TRUE, 1, 
                           FALSE, 1)
  }
  
  #drop infinite and nan propensities upfront
  dat <- filter(dat, !is.na(prop_score), is.finite(prop_score))
  
  #First reorganize the data by (t, k) and limit to interesting things
  #Keep track of number dropped for reporting
  q_thresh = quantile(dat$prop_score, q_prop_thresh)
  dropped_stats <- dat %>% mutate(to_drop = prop_score > q_thresh, 
                 numTests_cf = prob_counterfactual_test * estNumArrivals) %>%
    group_by(to_drop) %>%
    summarise( estNumArrivals = sum(estNumArrivals), 
               numPropScores = n(), 
               numTests_eva = sum(numTestsPerformed),
               numTests_counterfactual = sum(numTests_cf)) %>%
    mutate(test_scaling = numTests_eva / numTests_counterfactual)
  
  dat <- filter(dat, prop_score <= q_thresh)
  print(q_thresh)


  ## Compute the means!
  dat <- mutate(dat, is_peak = date_entry <= END_OF_SEASON)
  eva_means <- dat %>% group_by(is_peak) %>%
    summarise( eva_mean  = sum(numPositivesFound))

  #compute the policy means
  #first compute the infections time-series for plotting purposes
  infections <- dat %>%
    mutate(perf = numPositivesFound * prop_score) %>%
    group_by(date_entry) %>%
    summarise(Eva = sum(numPositivesFound), 
              cf_policy = sum(perf)) %>% 
    mutate(is_peak = date_entry <= END_OF_SEASON)  
  
  #use the scalings to inflate these
  infections <- left_join(infections, test_scalings, by="is_peak") %>%
    mutate(cf_policy = cf_policy * scaling) %>%
    select(-is_peak, -scaling)

  

  #now compute the averages over the seasons with ratios. 
  policy_means <- dat %>%
    mutate(perf = numPositivesFound * prop_score) %>%
    group_by(is_peak) %>%
    summarise( eva_mean = sum(numPositivesFound),
               cf_mean = sum(perf)) %>%
    mutate(ratio = cf_mean / eva_mean)

  policy_means <- left_join(policy_means, test_scalings, by="is_peak") %>%
    mutate(scaled_ratio = ratio * scaling)
  
  #### Now start computing variance stuff  
  dat.tk <- dat %>% 
    select(date_entry, eb_type, point_entry, numTestsPerformed, 
           numPositivesFound, prop_score, 
    ) %>%
    group_by(date_entry, eb_type) %>%
    summarise(numTestsPerformed = sum(numTestsPerformed), 
              numPositivesFound = sum(numPositivesFound), 
              prop_score = mean(prop_score), #these are all same, so mean is trivial
    ) %>%
    mutate(r_naive = if_else(numTestsPerformed > 0, numPositivesFound / numTestsPerformed, 0),
           is_peak = date_entry <= END_OF_SEASON)
  
  #compute the second component.
  dat.t_1 <- 
    dat.tk %>% 
    mutate(sec_term = prop_score^2 * numTestsPerformed * r_naive * (1 - r_naive))  %>%
    group_by(date_entry) %>%
    summarise(sec_term = sum(sec_term))
  
  #Now reorganize the data by (t, p) and limit to interesting things to compute 
  #first term
  dat.tkp <- dat %>% 
    select(date_entry, eb_type, point_entry, 
           estNumArrivals, numTestsPerformed, numTestsAlloc, numPositivesFound, 
           metric
    )
  
  #First marry on the r(t,k) estimates.  Notice these are at grouping (t,k) not (t,p)  
  dat.tkp <- left_join(dat.tkp, select(dat.tk, date_entry, eb_type, r_naive), 
                       by=c("date_entry", "eb_type"))
  
  dat.tp <- dat.tkp %>%
    mutate( num = r_naive * estNumArrivals * metric, 
            denom = estNumArrivals * metric) %>%
    group_by(date_entry, point_entry) %>%
    summarise( num = sum(num), 
               denom = sum(denom), 
               numTestsPerformed = sum(numTestsPerformed), 
               numTestsAlloc = sum(numTestsAlloc), 
               Ntot = sum(numTestsPerformed), 
    ) %>%
    mutate(w = if_else(denom > 0, num/denom, 0), 
           qns = if_else(numTestsAlloc > 0, 1- numTestsPerformed/numTestsAlloc, 0)
    )
  rm(dat.tkp)
  
  #roll it down to t level.
  dat.t_2 <- dat.tp %>%
    mutate(first_term = w^2 * numTestsAlloc * qns * (1 - qns)) %>%
    group_by(date_entry) %>%
    summarise(first_term = sum(first_term))
  
  dat.t <- left_join(dat.t_1, dat.t_2)
  rm(dat.t_1, dat.t_2)
  
  vs <- dat.t %>%
    mutate( mod_row = row_number() %% DELTA_AC, #computes   modulo arithmetic
            v_yt = first_term + sec_term, 
            is_peak = date_entry <= END_OF_SEASON ) %>%
    group_by(is_peak, mod_row) %>%
    summarise(var = sum(v_yt)) %>%
    mutate(std = sqrt(var))
  
  cf_std <- vs %>% group_by(is_peak) %>%
     summarise( std = sqrt(sum(outer(std,std))))
  
  #marry on the mean stuff for scaling.  
  cf_std <- left_join(cf_std, select(policy_means, -cf_mean)) %>%
    mutate(ratio = std/eva_mean) %>%
    select(-eva_mean) %>%
    mutate(scaled_ratio = ratio * scaling)

  return(list(infections=infections, 
              dropped_stats = dropped_stats,
              policy_means = policy_means, 
              cf_std = cf_std) )
}

#For the random policy
dat <- read_csv("../OPE_outputs/ope_dat_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv") %>%
  mutate(prob_counterfactual_test = prob_rand_test, 
         metric = 1)

compute_var_auc(dat)
rm(dat)

#load up the scaling data from Hamsa.
test_scalings <- read_tsv("../OtherData/testScalings.csv") %>%
  rename(peak=`Peak Scaling`, 
         off_peak = `OffPeak Scaling`) %>%
  pivot_longer(-c("Metric")) %>% 
  mutate(is_peak = name=="peak") %>% 
  select(-name) %>% 
  rename(scaling=value)


##Positivity data
dat <- read_csv("../OPE_outputs/ope_dat_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv")
dat.metric_prob <- read_csv("../OtherData/prob_public_pos_test.csv") %>%
  rename(date_entry=date, 
         prob_counterfactual_test = prob_public_pos_test)
dat.metric_metric <- read_csv("../OtherData/pos_metric.csv") %>%
  rename(date_entry=date, 
         metric = pos_metric)

dat.metric <- left_join(dat.metric_prob, dat.metric_metric)
rm(dat.metric_prob,dat.metric_metric)

dat.metric <- left_join(dat, dat.metric, by=c("date_entry", "country"="country"))

dat.metric <- dat.metric %>% 
  mutate(prop_score = prob_counterfactual_test / prob_gittins_test)

t.test_scalings <- filter(test_scalings, Metric=="Positivity")  %>%
  select(-Metric)

out <- compute_var_auc(dat.metric, test_scalings = t.test_scalings)
infections.pos <- out$infections %>% rename(Positivity=cf_policy)
out$policy_means
out$cf_std

rm(dat, dat.metric, out, t.test_scalings)

#Same game for Deaths
dat <- read_csv("../OPE_outputs/ope_dat_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv")
dat.metric_prob <- read_csv("../OtherData/prob_public_deaths_test.csv") %>%
  rename(date_entry=date, 
         prob_counterfactual_test = prob_public_deaths_test)
dat.metric_metric <- read_csv("../OtherData/deaths_metric.csv") %>%
  rename(date_entry=date, 
         metric = deaths_metric)

dat.metric <- left_join(dat.metric_prob, dat.metric_metric)
rm(dat.metric_prob,dat.metric_metric)

dat.metric <- left_join(dat, dat.metric, by=c("date_entry", "country"="country"))

dat.metric <- dat.metric %>% 
  mutate(prop_score = prob_counterfactual_test / prob_gittins_test)

t.test_scalings <- filter(test_scalings, Metric=="Deaths")  %>%
  select(-Metric)

out <- compute_var_auc(dat.metric, test_scalings = t.test_scalings)
infections.deaths<- out$infections %>% rename(Deaths=cf_policy)
out$policy_means
out$cf_std
rm(dat, dat.metric, out, t.test_scalings)

#once more for cases
dat <- read_csv("../OPE_outputs/ope_dat_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv")
dat.metric_prob <- read_csv("../OtherData/prob_public_cases_test.csv") %>%
  rename(date_entry=date, 
         prob_counterfactual_test = prob_public_cases_test)
dat.metric_metric <- read_csv("../OtherData/cases_metric.csv") %>%
  rename(date_entry=date, 
         metric = cases_metric)

dat.metric <- left_join(dat.metric_prob, dat.metric_metric)
rm(dat.metric_prob,dat.metric_metric)

dat.metric <- left_join(dat, dat.metric, by=c("date_entry", "country"="country"))

dat.metric <- dat.metric %>% 
  mutate(prop_score = prob_counterfactual_test / prob_gittins_test)

t.test_scalings <- filter(test_scalings, Metric=="Cases")  %>%
  select(-Metric)

out <- compute_var_auc(dat.metric, test_scalings = t.test_scalings)
infections.cases<- out$infections %>% rename(Cases=cf_policy)
out$policy_means
out$cf_std
rm(dat, dat.metric, out, t.test_scalings)

### 
# Now put it all together to make the pretty plot
####
infections <- rbind(pivot_longer(infections.cases, -date_entry, names_to="Metric"), 
                    pivot_longer(select(infections.deaths, -Eva), -date_entry, names_to="Metric"), 
                    pivot_longer(select(infections.pos, -Eva), -date_entry, names_to="Metric") 
              )
infections<- infections %>% pivot_wider(names_from = Metric)
write_csv(infections, "../OtherData/dailyTimeSeries_PublicPolicies.csv")

