######
# Generates historical time series of EB estimates
# Recall: In bandit, at time t we make a prediction for time t+1 using data [t-16, t-2].  Equivalently,
#         the time t prediction uses data [t-17, t-3].
#         eb_prev below uses data [t- prev_window, t + prev_window] for form estimate at time t
######
library(tidyverse)
library(lubridate)
source("mse_helpers.R")

#MIN_TEST  - Only countries with at least thsi many tests are used in the prior fiting step
#prev_window - Data from time [t-prev_window, t+prev_window] used for estimating prev at time t
#USE_MM - should moment-matching or MSE minimization be used to fit priors
run_file <- function(MIN_TEST, prev_window, USE_MM){

  ##CONSTANTS
  HIST_START_DATE = ymd("2020-08-02")  #System went live 5 Aug, but we are confident about the logging from the 2nd
  hist_end_date   = ymd("2020-11-03")  #Drops last 48 hours for testing delay.  We will only really go to 1 Nov.
  
  #Load and Clean Data
  plf_path <- str_c("../sample_data.csv")
  plf_data <- read_csv(
                       plf_path, 
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
  
  plf_data <- clean_hist_plf_data(plf_data, date_entry, "temp_log.log", 
                                  use_city_types=FALSE)
  plf_data <- filter(plf_data, date_entry <= hist_end_date)
  
  #Type labeling above done w.r.t current White List  
  #This is incorrect for USA which entered the whitelist on 17 Aug 2020.  Fix manually.
  plf_data <- 
    plf_data %>% 
    mutate(isCtryFlagged = ifelse(country == "US" & date_entry < ymd("2020-08-17"), 
                                  TRUE, isCtryFlagged))
  
  #assign colors for ease later
  plf_data <- plf_data %>%
    mutate(color = if_else(isCtryFlagged, "black", 
                           if_else(isCtryGrey, "grey", "white"))) 
  
  #####
  ##The major work.
  #Loop through history and compute the estimates for each day
  #####
  ix_today_dt = HIST_START_DATE
  tot_estimates <- NULL
  prior_estimates <-NULL
  
  while(ix_today_dt <= hist_end_date){
    #Create a temporary dataset with data from [t - prev_window, t + prev_window]
    ix_plf_data <- filter(plf_data, date_entry >= ix_today_dt - prev_window, 
                          date_entry <= ix_today_dt + prev_window)    
    hist_data <- 
      ix_plf_data %>% 
      group_by(eb_type, country, isCtryFlagged, isCtryGrey, color) %>%
      summarise( 
        num_arrivals = n(), 
        num_tested = sum(sent_for_test), 
        num_pos = sum(test_result == "positive", na.rm=TRUE), 
        num_inconclusive = sum(!is.na(test_result) & 
                               !test_result %in% c("positive", "negative"))
      )
    
    #Ensure historical data has entries for everyone on today's arrival list "pass_manifest"
    ix_pass_manifest <- 
      filter(plf_data, date_entry == ix_today_dt) %>%
      select(eb_type, isCtryFlagged, isCtryGrey, color)
    
    hist_data <- 
      full_join(
        hist_data, 
        unique(ix_pass_manifest), 
        by=c("eb_type", "isCtryFlagged", "isCtryGrey", "color") 
      ) %>%
      mutate(
        num_arrivals = ifelse(is.na(num_arrivals), 0, num_arrivals), 
        num_pos      = ifelse(is.na(num_pos), 0, num_pos), 
        num_tested   = ifelse(is.na(num_tested), 0, num_tested), 
        num_inconclusive = ifelse(is.na(num_inconclusive), 0, num_inconclusive)
      )

    #Debugging check only    
    if ( any(is.na(hist_data$isCtryFlagged)) | 
         any(is.na(hist_data$isCtryGrey)) | 
         any(is.na(hist_data$color))
    ){
      print("Problematic date ", ix_today_dt)
      break
    }
    
    ######
    #Now fit estimates.
    #Fit the EB models in simplest way possible for now.
    hist_data <- mutate(hist_data, prev = num_pos/num_tested) %>% ungroup() 

    if (USE_MM){
      #fit separate priors for each color by moment matching
      mom.by.color <- 
        hist_data %>% 
        filter(num_tested >= MIN_TEST) %>%
        group_by(color) %>%
        summarise(mom1 = mean(prev, na.rm = TRUE),
                  mom2 = mean(prev * (num_pos - 1)/(num_tested - 1), na.rm = TRUE))
      
      #MM may fail.  Correct by hand with appropriate default values
      #These defaults were chosen by "smoothing" adjacent priors where the fittin succeeded
      if (ix_today_dt <= ymd("2020-08-15")){
        mom_defaults <- 
          tribble(
            ~color, ~def_mom1, ~def_mom2, 
            "white",  0.002320979,   0.000011218, 
            "black",  0.00781,       0.000148, 
            "grey",   0.00171,       0.00000613
          )
      } else {
        mom_defaults <- 
          tribble(
            ~color, ~def_mom1, ~def_mom2, 
            "white",  0.002320979,   0.000011218, 
            "black",  0.00781,       0.000148, 
            "grey",   0.00243,       0.0000102
          ) }
      
      #Update moments for places where moments were not valid
      mom.by.color <- 
        left_join(mom.by.color, mom_defaults) %>%
        mutate(
          not_viable = non_viable_moments(mom1, mom2),
          mom1 = if_else(not_viable, def_mom1, mom1), 
          mom2 = if_else(not_viable, def_mom2, mom2)
        ) %>% 
        select(color, mom1, mom2)
      
      hist_data <- 
        left_join(hist_data, mom.by.color, by="color") %>%
        fit_eb_MM(mom1, mom2, MM)
      
    } else {   #Not using MM to fit priors. Using the MSE minimization instead.
      #Fit the mse priors
      hist_data <- mutate(hist_data, unit_weight = 1) 
      prior_by_color <- 
        map_df(
          c("white", "grey", "black"), 
          ~filter(hist_data, color ==., num_tested > MIN_TEST) %>%
            fit_p0_S0_mse(num_pos, num_tested, unit_weight) 
        )
      prior_by_color$color <- c("white", "grey", "black")
      
      prior_by_color <- 
        prior_by_color %>% 
        mutate(alpha = p0 * S0, beta = (1 - p0) * S0) %>%
        nest(MM = c(alpha, beta, p0, S0)) %>%
        select(-mse)
      
      hist_data <- left_join(hist_data, prior_by_color, "color")  
    }

    #Now equipped with daily priors, compute posteriors for each type.    
    hist_data <- 
      hist_data %>%
      add_eb_preds(eb_prev, MM, num_pos, num_tested) %>%
      unnest(MM) %>%
      mutate(h = alpha.post / (alpha.post + beta.post - 1), 
             hm = (alpha.post - 1) / (alpha.post + beta.post - 1))
    
    # Generate the current_estimates output
    curr_estimates <- 
      hist_data %>% 
      select(eb_type, country, isCtryFlagged, isCtryGrey, color, 
             eb_prev, prev, h, hm, num_pos, num_tested, num_arrivals, 
             alpha, beta, alpha.post, beta.post, p0, S0) %>%
      mutate(
        low = qbeta(.05, alpha.post, beta.post), 
        up  = qbeta(.95, alpha.post, beta.post), 
        date = ix_today_dt
      )
    
    #Add it to the output!
    tot_estimates = rbind(tot_estimates, curr_estimates)
    
    #increment for the while counter
    ix_today_dt = ix_today_dt + 1
  }  #loop to the next day
  
  #write it out.  Encode parameter settings in the name.  
  write_csv(
    tot_estimates,
    str_c("../OPE_outputs/hist_eb_timeseries_",
          USE_MM, "_Window_",if_else(is.na(prev_window), "NA", as.character(prev_window)), "_MinTest_", MIN_TEST, "_UnWeighted.csv")
  )
  
}

#These were the parameters used for the analysis in paper.  
#This output later smoothed -- see file smoothingEBPriors.r
run_file(30, prev_window = 3, USE_MM = TRUE)



