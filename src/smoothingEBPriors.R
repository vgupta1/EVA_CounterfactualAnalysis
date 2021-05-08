####  Smoothing EB Estimates
#  Reads in a generated EB File and does some smoothing via cubic splines
library(tidyverse)
library(lubridate)
library(broom)
source("helpers.R")

#USE_THRESH:
# Flag whether to sue absolute thresholding (THRESH) or quantile thresholding (Q_THRESH)
#THRESH/Q_THRESH:  
#  When smoothing priors, considers any prior where strength in prior exceeds a threshold.  
#  if USE_THRESH= TRUE, threshold is THRESH.  Else, quantile(prior_strengths, Q_THRESH)
#s_input_file: 
#  file output by the genHistoricalEB_TS.R script.
#  In production run, this is:  "../OPE_outputs/hist_eb_timeseries_TRUE_Window_3_MinTest_30_UnWeighted.csv"
#Default values indicate what was used in paper.
smoothEBFits <- function(s_input_file, USE_THRESH = TRUE, THRESH = 2001, Q_THRESH=.9){
  dat.all <- read_csv(s_input_file) 
  
  dat.priors <- 
    dat.all %>%
    group_by(date, color) %>% 
    slice_head() %>% 
    select(-eb_type, -country, -isCtryFlagged) %>%
    ungroup()
  
  #Identify which priors were the result of a default value
  #Currently done in a lazy/non-Robust way:
  #These tables copied from genHistEbtimeSeries and MUST match those tables for safety
  mom_defaults1 <- 
    tribble(
      ~color, ~def_mom1, ~def_mom2, 
      "white",  0.002320979,   0.000011218, 
      "black",  0.00781,       0.000148, 
      "grey",   0.00171,       0.00000613
    )
  
  mom_defaults1 <- 
    fit_eb_MM(mom_defaults1, def_mom1, def_mom2, MM) %>%
    unnest(MM) %>% 
    select(color, p0, S0) %>%
    rename(p0_def = "p0", S0_def = "S0")
  
  mom_defaults2 <- 
    tribble(
      ~color, ~def_mom1, ~def_mom2, 
      "white",  0.002320979,   0.000011218, 
      "black",  0.00781,       0.000148, 
      "grey",   0.00243,       0.0000102
    )
  
  mom_defaults2 <- 
    fit_eb_MM(mom_defaults2, def_mom1, def_mom2, MM) %>%
    unnest(MM) %>% 
    select(color, p0, S0) %>%
    rename(p0_def = "p0", S0_def = "S0")
  
  #This is dangerous with floating point arithmetic, but stable enough.
  dat.priors <- 
    left_join(dat.priors, mom_defaults1, by="color") %>%
    mutate(isDefault = p0 == p0_def) %>% 
    select(-contains("_def"))
  
  dat.priors <- 
    left_join(dat.priors, mom_defaults2, by="color") %>%
    mutate(isDefault = isDefault | (p0 == p0_def)) %>% 
    select(-contains("_def"))
  
  #iterate explicitly
  colors = c("white", "black", "grey")
  out <- NULL
  for (c in colors){
    dat.color <- filter(dat.priors, color == c)
    
    #smooth first moment based one everyone except default values  
    tdat.fit <-  filter(dat.color, !isDefault)
    tfit_p0 <- with(tdat.fit, smooth.spline(date, p0, cv = TRUE))
    
    #Add on the fits, truncated to [0,1]
    dat.color$p0_smooth <- 
      pmin(pmax(0, predict(tfit_p0, as.numeric(dat.color$date))$y), 1)
    
    #smooth second moment based on everyone except default and extreme values  
    tdat.fit <- filter(dat.color, !isDefault)
    thresh <- if_else(USE_THRESH, THRESH, quantile(tdat.fit$S0, Q_THRESH))
    tdat.fit <- filter(tdat.fit, S0 < thresh)
    tfit_S0 <- with(tdat.fit, smooth.spline(date, S0, cv = TRUE))
    
    #add on fits, avoiding zeros for stability.
    dat.color$S0_smooth <- pmax(0.0001, 
                                predict(tfit_S0, as.numeric(dat.color$date))$y)
    
    out <- rbind(out, dat.color)
  }
  dat.priors <- out
  rm(out)
  
  ##use these new priors to update the original posteriors
  dat.all <- left_join(dat.all, 
                       select(dat.priors, date, color, p0_smooth, S0_smooth))
  
  #Remove any of the old estimators to avoid confusion
  #update the p0, S0 and remove cruft
  #add back alpha, beta, alpha.post, beta.post, eb_prev
  dat.all <- 
    dat.all %>% 
    mutate(eb_prev_old = eb_prev) %>%
    select(-eb_prev, -h, -hm, -alpha, -beta, -alpha.post, -beta.post, -low, -up) %>%
    mutate(p0 = p0_smooth, S0=S0_smooth) %>% 
    select(-p0_smooth, -S0_smooth) %>%
    mutate(alpha = p0 * S0,
           beta  = (1 - p0) * S0, 
           alpha.post = alpha + num_pos, 
           beta.post = beta + num_tested - num_pos, 
           eb_prev = alpha.post / (alpha.post + beta.post), 
           h = alpha.post / (alpha.post + beta.post - 1), 
           hm = (alpha.post - 1) / (alpha.post + beta.post - 1)
    )
  
  #construct the output file name 
  #Retains everything but "unweighted" and tack on parameter choices into string
  out_file <- str_remove(s_input_file, "UnWeighted.*")
  out_file <- str_c(out_file, "SmoothPrior_", USE_THRESH, "_", 
                    THRESH, "_", Q_THRESH, ".csv")
  
  dat.all %>% write_csv(out_file)
}

#This was the run used for the analysis in paper
smoothEBFits("../OPE_outputs/hist_eb_timeseries_TRUE_Window_3_MinTest_30_UnWeighted.csv")
