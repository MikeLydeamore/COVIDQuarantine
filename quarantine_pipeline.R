library(reticulate)
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)

use_condaenv("./python/envs")

source_python("./python/quarantine_pipeline.py")
source_python("./python/infectionspread.py")

# % of concurrent outbreaks
tps <- c(1.2, 1.3, 1.5, 1.6, 1.7, 1.8)
quarantine_durations <- c(7L, 14L)

outbreak_probs <- lapply(tps, function(tp) {
  print(paste0("Running TP=", tp))
  lapply(quarantine_durations, function(duration) {
    x <- concurrent_outbreak_probabilities(tp, breaches_over_horizon=py_none(), quarantine_type="HomeQuar", quarantine_duration=duration)  
    y <- data.table(prob=as.vector(x))[, concurrent_outbreaks := 0:(.N-1)][, tp := tp][, duration := duration]
    
    return (y)  
  }) %>% rbindlist()
  
}) %>% rbindlist()

time <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
saveRDS(outbreak_probs, file=paste0("../results/outbreak_probs_", time, ".rds"))

facet_lookup <- paste0(quarantine_durations, " Day Quarantine")
names(facet_lookup) <- quarantine_durations

p <- ggplot(outbreak_probs, aes(x=concurrent_outbreaks, y=prob, fill=as.factor(tp))) + facet_wrap(~duration, labeller = as_labeller(facet_lookup)) +
  geom_col(position="dodge", colour="black") +
  labs(x="Number of concurrent outbreaks", y="Probability", fill="Transmission Potential") +
  theme_classic() + cowplot::panel_border("black") +
  scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks=seq(0, max(outbreak_probs$concurrent_outbreaks)),
                   labels=seq(0, max(outbreak_probs$concurrent_outbreaks)))

cowplot::save_plot(paste0("../results/outbreak_probs_", time, ".png"), type="cairo", plot=p, ncol=2)


# Time to major outbreak
probability_of_major_outbreak <- function(max_iterations, prob_threshold, transmission_potential, quarantine_type, quarantine_duration) {
  
  filename <- paste0("./data/penetration/TUnvacc_", quarantine_type, "_", quarantine_duration, "d_traveller_breach_timeseries.csv")
  df <- fread(filename)
  
  df[, time_left_quarantine := t_incubation + t_post_incubation - days_infectious_community]
  
  prob <- c(0)
  initial_infected <- 0L
  while (tail(prob, n=1) < prob_threshold) {
    initial_infected <- initial_infected + 1L
    cat(sprintf("\r"))
    message(paste0("Previous prob: ", tail(prob, n=1),". Now simulating with ", initial_infected, " initially infected."), appendLF=F)
    x <- lapply(1:max_iterations, function(i) {
      run_outbreak(num_escaped=initial_infected, max_generations=Inf, generation_zero_times = py_none(), 
                   generation_zero_time_left_quarantine = c(sample(df$time_left_quarantine, replace=T, size=initial_infected), NA),
                   max_infections = 60L,
                   transmission_potential_mean = transmission_potential, dispersion = 0.1, 
                   in_isolation = py_none())
      
    })
    
    prob <- c(prob, sum(sapply(x, get_finalsize) > 20) / length(x))
    
  }
  message(paste0("\nCompleted with ", initial_infected, " initial infected. Probability of major outbreak: ", tail(prob, n=1), " (>", prob_threshold, ")"))
  
  ret <- data.table(initial_infected=1:length(prob), prob=prob)
  
  return (ret)
}

#To work out when a major outbreak will happen for a given # of infected, we just use the distribution of time between breaches.
timing_of_major_outbreak <- function(max_iterations, initial_infected, quarantine_type, quarantine_duration) {
  
  filename <- paste0("./data/penetration/TUnvacc_", quarantine_type, "_", quarantine_duration, "d_traveller_breach_timeseries.csv")
  
  df <- fread(filename)
  breach_time_distribution <- diff(df$time_discharged)
  
  breach_times <- replicate(n=max_iterations, { max(cumsum(sample(breach_time_distribution, size=initial_infected, replace=T))) } )
  
  return (breach_times)
}

ob_timings <- lapply(tps, function(tp) {
  print(paste0("Running TP=", tp))
  lapply(quarantine_durations, function(duration) {
    ob_prob <- probability_of_major_outbreak(1000, 0.95, tp, quarantine_type="HomeQuar", quarantine_duration=duration)[, quarantine_duration := duration]
    timings <- lapply(ob_prob$initial_infected, function(i) {
      timing <- timing_of_major_outbreak(10000, i, "HomeQuar", duration)
      
      ret <- data.table(lower95=quantile(timing, probs=0.025), median=median(timing), upper95=quantile(timing, probs=0.975), mean=mean(timing), 
                        initial_infected=i, quarantine_duration=duration)
      
      return (ret)
      
    }) %>% rbindlist()
    
    ob_prob_merged <- merge(ob_prob, timings, by=c('initial_infected', 'quarantine_duration'))
    ob_prob_merged[, c("timelower95", "timemedian", "timeupper95", "timemean") := list(today()+lower95, today()+median, today()+upper95, today()+mean)]
    
    ob_prob_merged[, transmission_potential := tp]
  }) %>% rbindlist()
})


ggplot(ob_prob_merged, aes(y=prob, colour=as.factor(transmission_potential))) + 
  geom_line(aes(x=timemedian)) + 
  geom_ribbon(aes(xmin=timelower95, xmax=timeupper95), alpha=0.1) +
  facet_wrap(~duration, labeller=as_labeller(facet_lookup)) +
  coord_cartesian(xlim=c(as.Date(NA), ymd("2022-12-31"))) +
  labs(x="Date", y="Prob. of major outbreak", colour="Transmission potential") +
  theme_classic() + cowplot::panel_border("black")
max_iterations <- 1000

