library(reticulate)
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)

use_condaenv("./python/envs")

source_python("./python/quarantine_pipeline.py")
source_python("./python/infectionspread.py")

# % of concurrent outbreaks
tps <- data.table(tp=c(1.2, 1.7, 2.4), gen0vac=c(0.8, 1.0, 1.272), gen0unvacc=c(1.8, 1.5, 2.4))
quarantine_durations <- c(7L, 14L)
vaccination_types <- c("Unvacc", "VaccAZ1", "VaccAZ2")
infected_proportions <- c(0.01, 0.001, 0.0001)
num_travellers <- c(3000, 6000, 12000)

outbreak_probs <- lapply(tps, function(tp) {
  print(paste0("Running TP=", tp))
  lapply(quarantine_durations, function(duration) {
    lapply(vaccination_types, function(vacc) {
      x <- concurrent_outbreak_probabilities(tp, breaches_over_horizon=py_none(), quarantine_type="HotelQuar", quarantine_duration=duration,
                                             travellers_vaccinated=vacc, max_iterations=1000L)  
      y <- data.table(prob=as.vector(x))[, concurrent_outbreaks := 0:(.N-1)][, tp := tp][, duration := duration][, vaccination_type := vacc]
      
      return (y)  
    }) %>% rbindlist()
    
  }) %>% rbindlist()
  
}) %>% rbindlist()

time <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
saveRDS(outbreak_probs, file=paste0("../results/outbreak_probs_", time, ".rds"))

facet_lookup <- c(paste0(quarantine_durations, " Day Quarantine"), "Unvaccinated", "AZ 1 dose", "AZ 2 dose")
names(facet_lookup) <- c(quarantine_durations, vaccination_types)

p <- ggplot(outbreak_probs[tp == 2.4], aes(x=concurrent_outbreaks, y=prob, fill=as.factor(tp))) + facet_grid(vaccination_type~duration, labeller = as_labeller(facet_lookup)) +
  geom_col(position="dodge", colour="black") +
  labs(x="Number of concurrent outbreaks", y="Probability") +
  guides(fill=F) +
  theme_classic() + cowplot::panel_border("black") +
  scale_y_continuous(expand=expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks=seq(0, max(outbreak_probs$concurrent_outbreaks)),
                   labels=seq(0, max(outbreak_probs$concurrent_outbreaks)))

cowplot::save_plot(paste0("./results/outbreak_probs_", time, ".png"), type="cairo", plot=p, ncol=2)




ob_probs <- lapply(tps$tp, function(tp) {
  lapply(quarantine_durations, function(duration) {
    lapply(vaccination_types, function(vacc) {
      message(paste0("Running TP=", tp, ", quarantine duration=", duration, ", vaccination type=", vacc, "."))
      p <- probability_of_major_outbreak(100000, tp, quarantine_type="HotelQuar", quarantine_duration=duration, travellers_vaccinated=vacc)
      
      data.table(ob_prob=p, transmission_potential=tp, quarantine_duration=duration, travellers_vaccinated=vacc)
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()

ob_timings <- lapply(tps$tp, function(tp) {
  ob_times <- lapply(quarantine_durations, function(duration) {
    lapply(vaccination_types, function(vacc) {
      lapply(infected_proportions, function(prop) {
        lapply(num_travellers, function(n) {
          ob_prob <- ob_probs[transmission_potential == tp & quarantine_duration == duration & travellers_vaccinated == vacc, ob_prob]
          message(paste0("Running TP=", tp, ", quarantine duration=", duration, ", vaccination type=", vacc, ". OB prob: ", ob_prob))
          breach_times_outbreak <- timing_of_major_outbreak(1e6, ob_prob=ob_prob, 
                                                            quarantine_type="HotelQuar", quarantine_duration=duration, 
                                                            travellers_vaccinated=vacc,
                                                            infected_proportion=prop,
                                                            num_arrivals=n)
          
          ret <- data.table(breach_times=breach_times_outbreak, 
                            transmission_potential=tp, 
                            quarantine_duration=duration, 
                            travellers_vaccinated=vacc,
                            infected_proportion=prop,
                            num_arrivals=n)
          
          return (ret)
        }) %>% rbindlist()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()

ob_probs_hotel <- lapply(tps$tp, function(tp) {
  message(paste0("Running TP=", tp, "."))
  p <- probability_of_major_outbreak(100000, tp, quarantine_type="HotelQuar", quarantine_duration=14L, travellers_vaccinated="Unvacc", type="worker")
      
  data.table(ob_prob=p, transmission_potential=tp, quarantine_duration=14L, travellers_vaccinated="Unvacc")
}) %>% rbindlist()
  

ob_timings_combined6 <- lapply(quarantine_durations, function(duration) {
  lapply(vaccination_types, function(vacc) {
    lapply(infected_proportions, function(prop) {
      lapply(num_travellers, function(n) {
        lapply(tps$tp, function(tp) {
          ob_prob <- c(ob_probs[transmission_potential == tp & quarantine_duration == duration & travellers_vaccinated == vacc, ob_prob],
                       ob_probs_hotel[transmission_potential == tp & quarantine_duration == 14L & travellers_vaccinated == "Unvacc", ob_prob])
          
          breach_times_outbreak <- timing_of_major_outbreak(1e6, ob_prob=ob_prob, 
                                                            quarantine_type="HotelQuar", quarantine_duration=duration, 
                                                            travellers_vaccinated=vacc,
                                                            infected_proportion=prop,
                                                            num_arrivals=n)
          
          ret <- data.table(breach_times=breach_times_outbreak,
                            transmission_potential=tp, 
                            quarantine_duration=duration, 
                            travellers_vaccinated=vacc,
                            infected_proportion=prop,
                            num_arrivals=n)
          
          return (ret)
        }) %>% rbindlist()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
})
  
summarised_timings_combined <- ob_timings_combined2[, .(chance50=quantile(ecdf(breach_times), probs=0.5), chance95=quantile(ecdf(breach_times), probs=0.95)), 
                                          by=list(quarantine_duration, transmission_potential, travellers_vaccinated, infected_proportion, num_arrivals)]
summarised_timings_combined

time <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
saveRDS(ob_timings, file=paste0("./results/outbreak_timings_", time, ".rds"))

ggplot(ob_timings, aes(x=breach_times, colour=as.factor(transmission_potential), fill=as.factor(transmission_potential))) + 
  stat_ecdf(pad=F) + 
  #geom_ribbon(aes(xmin=timelower95, xmax=timeupper95), colour=NA, alpha=0.1) +
  facet_grid(travellers_vaccinated~quarantine_duration, labeller=as_labeller(facet_lookup)) +
  coord_cartesian(xlim=c(NA, 365*2)) +
  labs(x="Days", y="Prob. of major outbreak", colour="Transmission potential", fill="Transmission potential") +
  theme_classic() + cowplot::panel_border("black")

cowplot::save_plot(filename=paste0("./results/outbreak_timings_", time, ".png"), plot=last_plot(), nrow=2, type="cairo")

summarised_timings <- ob_timings[, .(chance50=quantile(ecdf(breach_times), probs=0.5), chance95=quantile(ecdf(breach_times), probs=0.95)), 
                                 by=list(quarantine_duration, transmission_potential, travellers_vaccinated, infected_proportion, num_arrivals)]


summarised_timings[transmission_potential == 2.4 & travellers_vaccinated == "VaccAZ2" & quarantine_duration == 7,]
summarised_timings[infected_proportion == 0.01 & travellers_vaccinated == "VaccAZ2" & quarantine_duration == 7 & num_arrivals == 3000,]

#HomeQuarantine:
scenarios <- list.files("data/penetration", pattern="Run_AZ1*", full.names=T)
home_ob_probs4 <- lapply(2.4, function(tp) {
  pb <- txtProgressBar(min=0, max=length(scenarios), style = 3)
  ob_probs <- lapply(1:length(scenarios), function(i) {
    setTxtProgressBar(pb, value=i)
    fname <- paste0(scenarios[i], "/Traveller_breach.csv")
    
    ob_prob <- probability_of_major_outbreak(quarantine_type = "HomeQuar", quarantine_duration = 14L, travellers_vaccinated = NA,
                                             transmission_potential = tp, filename = fname, max_iterations = 500000)
    
    return (data.table(scenario=scenarios[i], ob_prob=ob_prob, transmission_potential=tp))
  }) %>% rbindlist()
  close(pb)
  
  return (ob_probs)
}) %>% rbindlist()

home_ob_timings <- lapply(1:length(scenarios), function(i) {
  ob_prob <- home_ob_probs4[transmission_potential == tp & scenario == scenarios[i], ob_prob]
  fname <- paste0(scenarios[i], "/Traveller_breach.csv")
  breach_times_outbreak <- timing_of_major_outbreak(1e6, ob_prob=ob_prob, 
                                                    quarantine_type="HomeQuar", quarantine_duration=14, 
                                                    travellers_vaccinated=vacc,
                                                    infected_proportion=0.01,
                                                    num_arrivals=3000,
                                                    filename=fname)
  
  
  ret <- data.table(breach_times=breach_times_outbreak, 
                    transmission_potential=tp, 
                    num_arrivals=3000,
                    scenario=scenarios[i])
  
  return (ret)
}) %>% rbindlist()

summarised_timings_home <- home_ob_timings[, .(chance50=quantile(ecdf(breach_times), probs=0.5), chance95=quantile(ecdf(breach_times), probs=0.95)), 
                                 by=scenario]
