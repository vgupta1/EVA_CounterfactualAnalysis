library(tidyverse)
source("helpers.R")

# Uses p0 and S0 specification to add eb predictions
# pos_name, total_name indicate names of the input data
# pred_name is where the posterior expectation goes
# param_name is a tbl containing the prior alpha, beta and posterior alpha, beta
add_eb_preds_by_param <- function(df, p0, S0, 
                                  pred_name, param_name, pos_name, total_name){
  #convert the p0, S0 spec to an alpha, beta Spec
  df %>% ungroup() %>%
    mutate(
      alpha      = S0 * p0, 
      beta       = S0 * (1 - p0),
      alpha.post = alpha + {{pos_name}}, 
      beta.post  = beta + {{total_name}} - {{pos_name}}, 
      "{{pred_name}}" := alpha.post / (alpha.post + beta.post) 
    ) %>%
    nest("{{param_name}}" := c(alpha, beta, alpha.post, beta.post))
}

#Computes the MSE for a LOO type estimator
binom_mse_by_param <- function(df, p0, S0, pos_name, total_name, weight_name){
  #call the above to add a column for the pred (named h)
  #Need to subtract one to make the LOO function work out properly
  t <- df %>% mutate(tot_m = {{total_name}} - 1) %>%
    add_eb_preds_by_param(p0, S0, h, tparam, {{pos_name}}, tot_m) %>%
    select(-tparam)

  # #call again to add a column for the prediction with one fewer positive (hm)
  t <- 
    mutate(t, Xm = {{pos_name}} - 1) %>%
    add_eb_preds_by_param(p0, S0, h_m, tparam, Xm, tot_m) %>%
    select(-tparam)
  
  #add weights as appropriate
  if (missing(weight_name)){
    t <- mutate(t, t_weight = 1)
  } else {
    t <- mutate(t, t_weight = {{weight_name}})
  }

  #now form the mse constituents
  out<- 
    t %>% 
    mutate( 
      tprev = {{pos_name}} / {{total_name}}, 
      mse = h^2 * (1- tprev) + h_m^2 * tprev - 2 * tprev * h_m, 
      mse = mse * t_weight
    ) %>% 
    filter({{total_name}} > 0) %>%
    summarise(mean(mse)) %>%
    as.numeric()

  return(out)
}

#Computes the MSE for a LOO type estimator
# prev_name is the column which contains raw prevalence for day (X/N)
# h and (resp. h_m) are the estimates for the day assuming X (resp. X-1) positives, N-1 Tests.
binom_mse <- function(df, prev_name, h_name, h_m_name, weight_name){
  #add weights as appropriate
  if (missing(weight_name)) {
    df <- mutate(df, t_weight = 1)
  } else {
    df <- mutate(df, t_weight = {{weight_name}})
  }

  #now form the mse constituents
  out <- 
    df %>% 
    mutate(
      tprev = {{prev_name}}, h = {{h_name}}, h_m = {{h_m_name}},
      mse = h^2 * (1 - tprev) + h_m^2 * tprev - 2 * tprev * h_m, 
      mse = mse * t_weight
    ) %>% 
    filter(!is.na(tprev)) %>%
    summarise(mean(mse)) %>%
    as.numeric()

  return(out)
}


#Compute closed-form mean given S0
#weight_name are for computing weighted mse
comp_p0 <- function(df, S0, pos_name, total_name, weight_name){
  if (missing(weight_name)){
    df <- mutate(df, t_weight = 1)
  } else {
    df <- mutate(df, t_weight = {{weight_name}})
  }

  out <- 0
  if (S0 <= 0){
    out <- 
      df %>% 
      filter({{total_name > 0}}) %>%
      mutate(tprev = {{pos_name}} / {{total_name}} ) %>%
      summarise( out = weighted.mean(tprev, t_weight)) %>%
      as.numeric()
  } else {  #Typical Case where S0 > 0
    out <- 
      df %>%
      filter({{total_name}} > 0) %>% 
      mutate( 
        tot_m = {{total_name}} - 1, 
        wght = t_weight * S0^2 / (S0 + tot_m)^2, 
        tprev  = {{pos_name}} / {{total_name}}
      ) %>%
      summarise(num = weighted.mean(tprev, wght)) %>%
      as.numeric() 
  }

  return(out) #default value
}

#Exhaustive search for best S0
#if S_grid is NULL, defaults to
# TRACE = TRUE causes function to return mse along entire s_grid as a df
fit_p0_S0_mse <- function(df, pos_name, tot_name, weight_name, 
                          S_grid=NULL, trace=FALSE){
  #If the S_grid is null, default to something
  if (is_null(S_grid)){
    S_grid <- c(1e-6, seq(1, 100, length.out=30), 
                seq(125, 500, length.out = 32), seq(600, 2000, length.out=20) )
  }  
  out <- 
    tibble(
    S0 = S_grid, 
    p0 = map_dbl(S_grid, 
                 ~comp_p0(df, ., {{pos_name}}, {{tot_name}}, {{weight_name}}))
    ) 

  out$mse <- 
    map2_dbl(out$S0, out$p0, 
             ~binom_mse_by_param(df, .y, .x, {{pos_name}}, {{tot_name}}, 
                                 {{weight_name}})
    )
  
  if (trace){
    return(out)
  }
  
  return(out[which.min(out$mse), ])
}


comp_mse_file <- function(file_name){
  dat <- 
    read_csv(file_name) %>%
    mutate(unit_weight = 1)

  out  <- binom_mse(dat, prev, h, hm, unit_weight)  
  out2 <- binom_mse(dat, prev, h, hm, num_arrivals)

  return(tibble(mse=out, wmse=out2))
}  