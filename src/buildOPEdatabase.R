#####
# Off-Policy + GCF Analysis of the bandit data
###
library(tidyverse)
library(lubridate)
source("mse_helpers.R")

#Input: 
#  Path to historical eb time series:  
#     e.g., "../OPE_Outputs/hist_eb_timeseries_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv"
#Output: 
#   High-level Summary statistics averaged over peak/non-peak time periods.
#   Writes a file with full ope_database
runOPE_Analysis <- function(str_path_to_eb_est){
  t <- str_split(str_path_to_eb_est, "hist_eb_timeseries", simplify = TRUE)
  str_path_output <- str_c(t[1,1], "ope_dat", t[1, 2]) 
  rm(t)
  
  log_file = "log_counterfactual.log"
  END_OF_TIME = ymd("2121-01-01")
  GREY_LAG = 9  #number of extra days random policy takes to figure identify grey-listing
  PROP_CUT_OFF = 34.7  #Chosen to drop at most 2.5% of arrivals.
  END_OF_SEASON = ymd("2020-10-01")  #cutoff between peak and off-peak travel
  
  BANDIT_LIVE_DATE <- ymd("2020-08-06") #system live on 5 aug, but in greek time-zone, hence 6.
  TEST_DATE = ymd("2020-11-01") #last date for which we have reliable test allocations and results from bandit
  
  ######
  # Loading data
  plf_path <- str_c("../sample_data.csv")
  plf_data <- read_csv(plf_path, 
                       col_types = cols( 
                         result_id = col_double(),
                         country = col_character(),
                         city = col_character(),
                         age = col_double(),
                         gender = col_character(),
                         date_entry = col_date(format = ""),
                         point_entry = col_character(),
                         created_at = col_datetime(format = ""),
                         to_test = col_character(),
                         test_result = col_character(),
                         created_at_lab_test = col_datetime(format = ""),
                         sent_for_test = col_double()
                       )
  )
  
  ##Types defined by country/color, not city. 
  #Also throw out anything past the last testing allocation date
  plf_data <- clean_hist_plf_data(plf_data, date_entry, log_file, use_city_types=FALSE)
  plf_data <- filter(plf_data, date_entry <= TEST_DATE, 
                               date_entry >= BANDIT_LIVE_DATE)   
  
  #Labeling above uses latest whitelist  
  #This labeling is incorrect for USA which entered the whitelist on 17 Aug 2020
  #correct by hand
  plf_data <- 
    plf_data %>% 
    mutate(isCtryFlagged = ifelse(country == "US" & date_entry < ymd("2020-08-17"), 
                                   TRUE, isCtryFlagged))
  
  #assign colors for ease later
  plf_data <- plf_data %>%
    mutate(color = if_else(isCtryFlagged, "black", 
                           if_else(isCtryGrey, "grey", "white") )
    ) 
  
  dat_test_results <- read_csv("../sample_test_results.csv")
  
  plf_data <-left_join(plf_data, dat_test_results, by = c("result_id" = "id"))
  rm(dat_test_results)
  
  ####
  #Some Filtering
  ##
  #@WARNING:  Occasionally, we see test_results for people bandit didn't ask for.
  #Thus, it's possible that numTestsPerformed > numTestsAlloc.
  #Omit such people and anyone "not_seen" since the algorithm didn't work on them
  plf_data <- filter(plf_data, !to_test == "not_seen",                   
                     test_was_allocated | is.na(test_result),  
  ) 

  #We model a "no-show" is someone Algorithm allocated a test but was not 
  # sent for test  (remember if flagged, then test_was_allocated = true)
  #Tacitly, this includes operational mishaps as "no-shows"
  #If no tests were allocated, we assume no-show rate is 0 below.
  plf_data <- 
    plf_data %>% 
    mutate(no_show =  ifelse(test_was_allocated, !sent_for_test, NA) ) 
  
  ######################
  #Melt down to summary statistics for convenience
  #compute no-show rate separately each port + type since at this level, 
  #bandit probabilities are constant
  ope_dat <- 
    plf_data %>% 
    group_by(date_entry, eb_type, country, color, point_entry) %>%
    summarise( numSchedArrivals  = n(), 
               numTestsPerformed = sum(!is.na(test_result)), 
               numPositivesFound = sum(test_result == "positive", na.rm = TRUE), 
               numTestsAlloc = sum(test_was_allocated),
               no_show_rate  = sum(no_show, na.rm = TRUE) / sum(!is.na(no_show)),
               no_show_rate  = ifelse(is.na(no_show_rate), 0., no_show_rate), #@WARNING:  there are NA's in the no-show rate, because we don't test anyone sometimes.  
               estNumArrivals = ifelse(numSchedArrivals > 0, 
                                       numSchedArrivals * (1 - no_show_rate), 0.)
               )
  
  #Read in eb_times series and marry on the estimates
  eb_ts <- 
    read_csv(str_path_to_eb_est) %>%
    select(-country, eb_prev, alpha.post, beta.post)

  ope_dat <- left_join(ope_dat, eb_ts,
                       by = c("date_entry"="date", "eb_type", "color"))
  
  #####
  # Labeling the Flux Periods for GreyListing
  # Assumes benchmark policy would grey-list a country GREY_LAG days after we did
  # Label all data points that occur in this window as "flux"
  ######
  grey_list_se <- 
    read_csv("../OtherData/grey_list_start_end.csv",
             col_types = cols(end_date = col_date()) ) %>%
    mutate(end_date = as_date( ifelse(is.na(end_date), END_OF_TIME, end_date) ))
  
  #Brute force-ish. 
  ope_dat <- 
    left_join(ope_dat, grey_list_se, by="country") %>%
    mutate(is_pt_flux = start_date <= date_entry & date_entry <= start_date + GREY_LAG,
           is_pt_flux = ifelse(is.na(is_pt_flux), FALSE, is_pt_flux), 
           eb_prev_flux = eb_prev,
           estNumArrivals_flux = estNumArrivals
    ) %>% ungroup()
  
  #Update the eb_prev_flux and estNumArrivals_flux for anyone in a flux period
  #Read in the the GCF predictions for these quantities
  #EB Preds can just be married on
  #For arrivals, split them up in the same fractions as the estNumArrivals were split.
  gcf_eb_preds <- 
    read_csv("../grey_eb_preds.csv") %>%
    select(-X1) %>% rename(date_entry = "date")

  #marry on the GCF_eb_prev
  ope_dat <- left_join(ope_dat, 
            select(gcf_eb_preds, country, date_entry, pred_prev)) %>% 
    mutate(eb_prev_flux = if_else(is_pt_flux, pred_prev, eb_prev_flux)) %>%
    select(-pred_prev)
  
  #First, compute total arrivals at from each country each day  
  #Marry on GCF predictions
  arriv.day <- 
    ope_dat %>% 
    group_by(date_entry, country, is_pt_flux) %>%
    summarise(totEstNumArrivals = sum(estNumArrivals))
  
  gcf_arrivals <- 
    read_csv("../grey_arrivs.csv") %>% 
    select(-X1) %>% rename(date_entry = "date") %>%
    mutate(totPredArrivals = pred_arr)

  arriv.day <- 
    left_join(arriv.day, 
              select(gcf_arrivals, country, date_entry, totPredArrivals))
  
  #marry it back on so we can compute fractions
  t <- left_join(ope_dat, arriv.day)
  t <- 
    t %>% 
    mutate(
      frac_by_port_type = if_else(totEstNumArrivals > 0, 
                                  estNumArrivals / totEstNumArrivals, 0), 
      pred_arrivals = frac_by_port_type * totPredArrivals
    )
  
  ope_dat <- 
    t %>% 
    mutate(estNumArrivals_flux = if_else(is_pt_flux, pred_arrivals, 
                                         estNumArrivals_flux)) %>%
    select(-totEstNumArrivals, -totPredArrivals, 
           -frac_by_port_type, -pred_arrivals)
  
  #Compute the randomization probabilities for the random policy
  #We are assuming that Y (infection) independent of Port given X (covariates), i.e, port arrivals not informative
  #Under this assumption, we actually compute expected randomization across ports
  #We get two randomizations, based on if we use flux arrivals or not
  dat.by_port <- 
    ope_dat %>% 
    group_by(date_entry, point_entry) %>% 
    summarise(numTestsPerformed = sum(numTestsPerformed), 
              estNumArrivals = sum(estNumArrivals),
              estNumArrivals_flux = sum(estNumArrivals_flux),
              prob_test_at_port = ifelse(estNumArrivals > 0, 
                                         numTestsPerformed/estNumArrivals, 0), ###  Prob(Person Tested at Port  under Random) 
              prob_test_at_port_flux = ifelse(estNumArrivals_flux > 0, 
                                              numTestsPerformed/estNumArrivals_flux, 0)
    )
  
  dat.by_port_x <- 
    ope_dat %>% 
    group_by(date_entry, point_entry, eb_type) %>%
    summarise(estNumArrivals_port_x = sum(estNumArrivals), ###Number of arrivals at a given port of a given type
              estNumArrivals_port_x_flux = sum(estNumArrivals_flux) )
  
  dat.by_x <- 
    ope_dat %>% 
    group_by(date_entry, eb_type) %>%
    summarise(estNumArrivals_x = sum(estNumArrivals), ### Number of arrivals of given type
              estNumArrivals_x_flux = sum(estNumArrivals_flux) )
  
  dat.by_port_x <- 
    left_join(dat.by_port_x, dat.by_x, by = c("date_entry", "eb_type")) %>% 
    mutate( 
      prob_port_given_x = estNumArrivals_port_x / estNumArrivals_x, #Prob ( Arrive at Port = 1 | eb_type )
      prob_port_given_x = ifelse(estNumArrivals_x == 0, 0, prob_port_given_x),
      prob_port_given_x_flux = estNumArrivals_port_x_flux / estNumArrivals_x_flux, 
      prob_port_given_x_flux = ifelse(estNumArrivals_x_flux == 0, 
                                      0, prob_port_given_x_flux)
    )
  
  dat.by_x <- 
    left_join(dat.by_port_x, dat.by_port, by = c("date_entry", "point_entry")) %>% 
    mutate( prod = prob_test_at_port * prob_port_given_x,                             #Prob ( Arbitrary Person Tested at Port = 1) * Prob (Port = 1 | German ) = 
            prod_flux = prob_test_at_port_flux * prob_port_given_x_flux ) %>%
    group_by(date_entry, eb_type) %>%
    summarise( tau1 = sum(prod), 
               tau1_flux = sum(prod_flux))
  
  ope_dat <- left_join(ope_dat, dat.by_x, by=c("date_entry", "eb_type")) %>%
    rename("prob_rand_test" = "tau1", "prob_rand_test_flux" = "tau1_flux")
  
  rm(dat.by_x, dat.by_port, dat.by_port_x)
  
  #Compute the randomization probabilities for our policy
  dat.by_x <- 
    ope_dat %>% 
    group_by(date_entry, eb_type) %>%
    summarise(
      numTestsPerformed = sum(numTestsPerformed), 
      estNumArrivals = sum(estNumArrivals),
      prob_gittins_test = if_else( estNumArrivals > 0, 
                                   numTestsPerformed/estNumArrivals, 0), #@warning we can have NA if no tests performed and no arrivals
    )
  
  ope_dat <- left_join(ope_dat, 
                       select(dat.by_x, date_entry, eb_type, prob_gittins_test), 
                       by=c("date_entry", "eb_type"))
  rm(dat.by_x)
  
  ####
  ##Form Estimates without the Grey-listing Counterfactuals
  #Direct method - Likely systematically biased
  ope_dat <- 
    ope_dat %>% 
    mutate(
      estNumPositive     = estNumArrivals * eb_prev,
      rand_caught_direct = estNumArrivals * prob_rand_test * eb_prev, 
      gitt_caught_direct = estNumArrivals * prob_gittins_test * eb_prev
    )
  
  #IPW method, and the doubly robust method - Unbiased by construction
  ope_dat <- ope_dat %>% 
    mutate( 
      prop_score = prob_rand_test / prob_gittins_test, 
      rand_caught_ipw = ifelse(numPositivesFound > 0, 
                               numPositivesFound * prop_score, 0),  
      gitt_caught_ipw = numPositivesFound,
      rand_caught_dr  = numPositivesFound - eb_prev * numTestsPerformed, 
      rand_caught_dr  = rand_caught_dr * prop_score,
      rand_caught_dr  = ifelse(is.nan(prop_score) | prop_score > PROP_CUT_OFF, 
                               0., rand_caught_dr), #omit unstable propensities
      rand_caught_dr  = rand_caught_direct + rand_caught_dr, 
      gitt_caught_dr  = numPositivesFound - eb_prev * numTestsPerformed, 
      gitt_caught_dr  = gitt_caught_dr + 
                          eb_prev * prob_gittins_test * estNumArrivals
    )
  
  #####
  # Repeat above calculations for the Grey-list counterfactuals
  #Random allocations uses the flux versions of eb_prev, estNumArrivals, and prob_rand_test
  #Direct method - possibly systematically biased
  ope_dat <- 
    ope_dat %>% 
    mutate(
      estNumPositive_GCF = estNumArrivals_flux * eb_prev_flux,
      GCF_caught_direct = estNumArrivals_flux * prob_rand_test_flux * eb_prev_flux
    )
  
  #IPW method, and the doubly robust method don't apply out of box for GCF
  #because arrivals are not unbiased for the new prevalence.  

  #Add on variance approximation information
  ope_dat <- 
    ope_dat %>% 
    mutate(qk = if_else(estNumArrivals > 0, numPositivesFound / estNumArrivals, 0),
           var_contrib = prop_score^2 * estNumArrivals * qk * (1 - qk))
  
  ##Dump this complete table of OPE output for summary level analysis for paper  
  write_csv(ope_dat, str_path_output)
    
  #Now computes some top-level summary stats for fun
  ###Focus on peak season and reasonable propensities
  out.sum <- 
    ope_dat %>% ungroup() %>%
    mutate(is_peak = date_entry <= END_OF_SEASON) %>%
    group_by(is_peak) %>%
    filter(prop_score <= PROP_CUT_OFF) %>%
    select(
      estNumPositive, prop_score, var_contrib,
      contains("caught"), contains("missed"),
      numPositivesFound, numTestsAlloc, numTestsPerformed,
      estNumArrivals, numSchedArrivals
    ) %>%
    summarise(across(-prop_score, .fns=sum) )

  # Reorganize to be human readable
  out.ratios <- 
    out.sum %>% 
    select(is_peak, var_contrib, contains("caught")) %>%
    pivot_longer(-c(is_peak, var_contrib),  
                 names_to = c("Policy", "Method"),
                 names_pattern = "(.*)_caught_(.*)",
                 values_to = "NumPeople"
    )
  
  out.ratios <- 
    left_join(
      out.ratios, 
      out.ratios %>% filter(Policy == "rand") %>% select(-Policy, -var_contrib), 
      by=c("Method", "is_peak")
    ) %>%
    mutate(
      Ratio_to_Rand = NumPeople.x / NumPeople.y, 
      NumPeople = NumPeople.x, 
      path = str_path_output, 
      std = if_else(Policy == "rand" & Method == "ipw", sqrt(var_contrib), 0.)
    ) %>%
    select(-NumPeople.x, -NumPeople.y) %>%
    select(-path, -var_contrib, path)

  #Do a back of the envelope calc to figure out how many sick people you prevented via grey-listing
  #home-rate = effective prevalence among people who stayed home as a result of grey list
  #Notice this only makes sense at the eb_type level.
  num_prevented <- 
    ope_dat %>% mutate(is_peak = date_entry <= END_OF_SEASON) %>%
    filter(is_pt_flux) %>%
    group_by(date_entry, eb_type, is_peak) %>%
    summarise(
      across(contains("NumArrivals"), sum), 
      across(contains("eb_prev"), mean),  #mean is a dummy function.  all entries identical anyway
    ) %>%
    mutate(
      home_prev = eb_prev_flux * estNumArrivals_flux, 
      home_prev = home_prev - eb_prev * estNumArrivals, 
      home_prev = home_prev / (estNumArrivals_flux - estNumArrivals),
      home_prev = pmin(pmax(0, home_prev), 1),
      num_prevented = pmax(estNumArrivals_flux - estNumArrivals, 0) * home_prev, 
      num_prevented = num_prevented - pmax(0, estNumArrivals - estNumArrivals_flux) * eb_prev
    ) %>%
    group_by(is_peak) %>%
    summarise(num_prevented = sum(num_prevented))

  return(list(ope_ratios=out.ratios, num_prevented=num_prevented, ope_summary = out.sum, ope_dat = ope_dat))    
}

#For by hand calculations, note that Grey-list coutnerfactual prediction reports
#Peak SD: 13.47438; Off-peak SD: 4.533486 for the 9 day period
# for the noise of the estimated GCF benefit 
t <- runOPE_Analysis("../Ope_outputs/hist_eb_timeseries_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv")

#Performance Ratios
perf <-
  t$ope_summary %>%
  select(is_peak, gitt_caught_ipw, rand_caught_ipw, var_contrib) %>%
  mutate(ratio = gitt_caught_ipw / rand_caught_ipw,
         sd_ratio = sqrt(var_contrib) / rand_caught_ipw)

#Inspect the adjustments for GCF
left_join(perf, t$num_prevented) %>%
  mutate(
    prevented_ratio = num_prevented / rand_caught_ipw,
    GCF_sd_ratio = if_else(is_peak, 13.47438, 4.533486) / rand_caught_ipw, ###This uses the 9day window!!!!
    GCF_ratio = (gitt_caught_ipw + num_prevented)/ rand_caught_ipw)
