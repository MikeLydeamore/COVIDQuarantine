# Time to major outbreak
probability_of_major_outbreak <- function(max_iterations, transmission_potential,
                                          quarantine_type, quarantine_duration, travellers_vaccinated, type,
                                          filename=NA,
                                          tp, tpgen0, tpgen0_vaccinated) {
  
  if (is.na(filename)) {
    filename <- paste0("./data/penetration/T", travellers_vaccinated, "_WUnvacc_", quarantine_type, "_", quarantine_duration, "d_", type, "_breach_timeseries.csv")
  }
  df <- fread(filename)
  
  if (type == "traveller") {
    if (quarantine_type == "HotelQuar")
      df[, time_left_quarantine := t_incubation + t_post_incubation - days_infectious_community][, time_entered_isolation := Inf]
    else
      df[, time_left_quarantine := t_incubation + t_post_incubation - exposure_days][, time_entered_isolation := Inf]
  }
  
  else if (type == "worker")
    df[, time_left_quarantine := 0][, time_entered_isolation := days_to_detection]
  
  
  outbreaks <- replicate(n = max_iterations, { 
    
    sampled_row <- df[sample(.N, 1), ]
    gen0tp <- ifelse(sampled_row$vaccinated == 1, tpgen0_vaccinated, tpgen0)
    
    if (type == "traveller") {
      run_outbreak(num_escaped=1L, max_generations=Inf, generation_zero_times = py_none(), 
                   generation_zero_time_left_quarantine = c(sampled_row$time_left_quarantine, NA),
                   max_infections = 10L,
                   transmission_potential_mean = transmission_potential, generation_zero_transmission_potential_mean = gen0tp,
                   dispersion = 0.1, 
                   in_isolation = py_none()) 
    } else {
      run_outbreak(num_escaped=1L, max_generations=Inf, generation_zero_times = py_none(), 
                   generation_zero_time_left_quarantine = c(sampled_row$time_left_quarantine, NA),
                   max_infections = 10L,
                   transmission_potential_mean = transmission_potential, generation_zero_transmission_potential_mean = gen0tp,
                   dispersion = 0.1, 
                   in_isolation = TRUE, generation_zero_time_entered_isolation = c(sampled_row$time_entered_isolation, NA)) 
    }
    
  }, simplify = FALSE )
  
  prob <- sum(sapply(outbreaks, get_finalsize) >= 5) / length(outbreaks)
  
  #message(paste0("\nCompleted with ", 1, " initial infected. Probability of major outbreak: ", tail(prob, n=1), " (>", prob_threshold, ")"))
  
  return (prob)
}

#To work out when a major outbreak will happen for a given # of infected, we just use the distribution of time between breaches.
timing_of_major_outbreak <- function(max_iterations, ob_prob, quarantine_type, quarantine_duration, 
                                     travellers_vaccinated, num_arrivals, infected_proportion, filename=NA,
                                     scaling_factor=NA) {
  if (is.na(filename)) {
    filename <- paste0("./data/penetration/T", travellers_vaccinated, "_WUnvacc_", quarantine_type, "_", quarantine_duration, "d_traveller_breach_timeseries.csv")
  }
  df <- fread(filename)
  
  if (quarantine_type == "HomeQuar") {
    breach_time_distribution <- diff(df$time_discharged)
    if (is.na(scaling_factor)) {
      scaling_factor <- num_arrivals / 100 * (infected_proportion / 0.01) #100 arrivals from Cam's outputs
      if (quarantine_duration == 7) {
        scaling_factor <- scaling_factor * 2
      }
    }

    breach_time_distribution <- breach_time_distribution / scaling_factor
    
    breach_times <- sample(breach_time_distribution, size=max_iterations, replace=T)
    is_outbreak <- rbinom(n=length(breach_times), size=1, prob=ob_prob)
    
    breach_times_outbreak <- cumsum(breach_times)[is_outbreak == 1] %>% diff()
    
    return (breach_times_outbreak)
  }
  
  else {
    
    fname_workers <- paste0("./data/penetration/T", travellers_vaccinated, "_WUnvacc_", quarantine_type, "_", quarantine_duration, "d_worker_breach_timeseries.csv")
    df_workers <- fread(fname_workers)
    df_travellers <- df
    
    df_combined <- rbind(df_travellers[, .(time_discharged)], df_workers[, .(time_discharged = time_removed-days_to_detection)])
    df_combined <- df_combined[order(time_discharged)]
    
    breach_time_distribution_travellers <- diff(df_travellers$time_discharged)
    breach_time_distribution_workers <- diff(df_workers[, time_removed - days_to_detection])
    
    breach_time_distribution <- diff(df_combined$time_discharged)
    scaling_factor <- num_arrivals / 100 * (infected_proportion / 0.01) #100 arrivals from Cam's outputs
    if (quarantine_duration == 7) {
      scaling_factor <- scaling_factor * 2
    }
    
    breach_time_distribution_travellers <- breach_time_distribution_travellers / scaling_factor
    breach_time_distribution_workers <- breach_time_distribution_workers / scaling_factor
    breach_time_distribution <- breach_time_distribution / scaling_factor
    
    prop_travellers <- df_travellers[,.N] / (df_workers[,.N] + df_travellers[, .N])
    num_travellers_sample <- round(prop_travellers * max_iterations)
    
    breach_times_travellers <- sample(breach_time_distribution_travellers, size=num_travellers_sample, replace=T)
    breach_times_workers <- sample(breach_time_distribution_workers, size=max_iterations-num_travellers_sample, replace=T)
    breach_times <- sample(breach_time_distribution, size=max_iterations, replace=T)
    
    is_outbreak_travellers <- rbinom(n=length(breach_times_travellers), size=1, prob=ob_prob[1])
    is_outbreak_workers <- rbinom(n=length(breach_times_workers), size=1, prob=ob_prob[2])
    
    is_outbreak <- c(rbinom(n=num_travellers_sample, size=1, prob=ob_prob[1]), rbinom(n=max_iterations-num_travellers_sample, size=1, prob=ob_prob[2]))
    is_outbreak <- is_outbreak[sample(1:length(is_outbreak), size=length(is_outbreak), replace=F)]
    
    #indexes <- sample(x = 1:length(breach_times_all), size = length(breach_times_all), replace = FALSE)
    
    #breach_times_all <- breach_times_all[indexes]
    #is_outbreak_all <- is_outbreak_all[indexes]
    
    breach_times_outbreak <- cumsum(breach_times)[is_outbreak == 1] %>% diff()
    
    ret <- list(breach_times_workers = cumsum(breach_times_workers)[is_outbreak_workers == 1] %>% diff(), 
                breach_times_travellers = cumsum(breach_times_travellers)[is_outbreak_travellers == 1] %>% diff(),
                prop_travellers = prop_travellers)
    
    return (breach_times_outbreak)
    
  }
  
  
}
