library(reticulate)
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)

use_condaenv("./python/envs")

source_python("./python/quarantine_pipeline.py")
source_python("./python/infectionspread.py")

source("quarantine_pipeline_functions.R")

R0s <- 5:10
VEs <- seq(0.5, 0.9, by=0.1)
coverages <- seq(0.1, 1.0, by=0.1)

vars <- expand.grid(R0=R0s, VE=VEs, coverage=coverages)

ob_probs <- mapply(R0=vars$R0, VE=vars$VE, coverage=vars$coverage, FUN = function(R0, VE, coverage) {
  
  message(paste0("Running R0: ", R0, ", VE: ", VE, ", Coverage: ", coverage))
  
  vestring <- as.character(VE) %>% gsub("\\.", "-", .)
  
  fname1 <- paste0("./paper/fixed/Run_R0_", R0, "-0_VE_", vestring, "_2021_08_26")
  fname2 <- paste0("./paper/run_2_2021_08_06/Run_R0_", R0, "-0_VE_", vestring, "_2021_08_06")
  
  tp <- get_tps(R0, prop_vac = coverage, efficacy_infection = VE, efficacy_transmission = VE)
  
  traveller_ob_prob <- probability_of_major_outbreak(max_iterations = 100000, 
                                  transmission_potential = tp$new_pop_tp, 
                                  tpgen0 = tp$new_tp_import_unvac,
                                  tpgen0_vaccinated = tp$new_tp_import_vac,
                                  type = "traveller",
                                  filename = paste0(fname1, "/Traveller_breach.csv"),
                                  quarantine_type = "HomeQuar",
                                  travellers_vaccinated = "Vacc",
                                  quarantine_duration = 14L)
  
   worker_ob_prob <- probability_of_major_outbreak(max_iterations = 100000, 
                                   transmission_potential = tp$new_pop_tp, 
                                   tpgen0 = tp$new_tp_import_unvac,
                                   tpgen0_vaccinated = tp$new_tp_import_vac,
                                   type = "worker",
                                   filename = paste0(fname1, "/Worker_breach.csv"),
                                   quarantine_type = "HomeQuar",
                                   travellers_vaccinated = "Vacc",
                                   quarantine_duration = 14L)
  
  data.table(R0, VE, coverage, traveller_ob_prob, worker_ob_prob)
  
}, SIMPLIFY = F)

ob_probs_extra <- mapply(R0=3, VE=0, coverage=vars$coverage, FUN = function(R0, VE, coverage) {
  
  message(paste0("Running R0: ", R0, ", VE: ", VE, ", Coverage: ", coverage))
  
  vestring <- "0-0"
  
  fname1 <- paste0("./paper/fixed/Run_R0_", R0, "-0_VE_", vestring, "_2021_08_26")
  
  tp <- get_tps(R0, prop_vac = coverage, efficacy_infection = VE, efficacy_transmission = VE)
  
  traveller_ob_prob <- probability_of_major_outbreak(max_iterations = 100000, 
                                                     transmission_potential = tp$new_pop_tp, 
                                                     tpgen0 = tp$new_tp_import_unvac,
                                                     tpgen0_vaccinated = tp$new_tp_import_vac,
                                                     type = "traveller",
                                                     filename = paste0(fname1, "/Traveller_breach.csv"),
                                                     quarantine_type = "HomeQuar",
                                                     travellers_vaccinated = "Vacc",
                                                     quarantine_duration = 14L)
  
  worker_ob_prob <- probability_of_major_outbreak(max_iterations = 100000, 
                                                  transmission_potential = tp$new_pop_tp, 
                                                  tpgen0 = tp$new_tp_import_unvac,
                                                  tpgen0_vaccinated = tp$new_tp_import_vac,
                                                  type = "worker",
                                                  filename = paste0(fname1, "/Worker_breach.csv"),
                                                  quarantine_type = "HomeQuar",
                                                  travellers_vaccinated = "Vacc",
                                                  quarantine_duration = 14L)
  
  data.table(R0, VE, coverage, traveller_ob_prob, worker_ob_prob)
  
}, SIMPLIFY = F) %>% rbindlist()

ob_probs <- rbindlist(ob_probs)
fwrite(ob_probs, file="paper/ob_probs_travellers.csv")

ggplot(ob_probs, aes(x=coverage, y=traveller_ob_prob, colour=as.factor(R0))) + facet_grid(.~VE) + geom_line()
ggplot(ob_probs, aes(x=coverage, y=worker_ob_prob, colour=as.factor(R0))) + facet_grid(.~VE) + geom_line()

ob_timings <- lapply(1:nrow(vars), function(i) {

  ob_prob <- c(ob_probs[R0 == vars$R0[i] & VE == vars$VE[i] & coverage == round(vars$coverage[i], digits=1), traveller_ob_prob],
             ob_probs[R0 == vars$R0[i] & VE == vars$VE[i] & coverage == round(vars$coverage[i], digits=1), worker_ob_prob])
  
  vestring <- as.character(vars$VE[i]) %>% gsub("\\.", "-", .)
  
  fname <- paste0("./paper/fixed/Run_R0_", vars$R0[i], "-0_VE_", vestring, "_2021_08_26")
  
  breach_times_outbreak <- timing_of_major_outbreak(1e6, ob_prob=ob_prob, 
                                                    quarantine_type="HotelQuar", quarantine_duration=14L, 
                                                    travellers_vaccinated="Vacc",
                                                    infected_proportion=0.1,
                                                    num_arrivals=3000,
                                                    filename=paste0(fname, "/Traveller_breach.csv"),
                                                    filename_workers=paste0(fname, "/Worker_breach.csv"))
  
  
  data.table(breach_times=breach_times_outbreak,
             R0=vars$R0[i], VE=vars$VE[i], coverage=vars$coverage[i])
  
})

ob_timings <- rbindlist(ob_timings)

ob_timings_extra <- lapply(1:nrow(vars), function(i) {
  
  ob_prob <- c(ob_probs[R0 == vars$R0[i] & VE == vars$VE[i] & coverage == vars$coverage[i], traveller_ob_prob],
               ob_probs[R0 == vars$R0[i] & VE == vars$VE[i] & coverage == vars$coverage[i], worker_ob_prob])
  print(ob_prob)
  vestring <- "0-0"
  
  fname <- paste0("./paper/fixed/Run_R0_", vars$R0[i], "-0_VE_", vestring, "_2021_08_26")
  
  breach_times_outbreak <- timing_of_major_outbreak(1e6, ob_prob=ob_prob, 
                                                    quarantine_type="HotelQuar", quarantine_duration=14L, 
                                                    travellers_vaccinated="Vacc",
                                                    infected_proportion=0.1,
                                                    num_arrivals=3000,
                                                    filename=paste0(fname, "/Traveller_breach.csv"),
                                                    filename_workers=paste0(fname, "/Worker_breach.csv"))
  
  
  data.table(breach_times=breach_times_outbreak,
             R0=vars$R0[i], VE=vars$VE[i], coverage=vars$coverage[i])
  
})

ob_timings <- rbind(ob_timings, rbindlist(ob_timings_extra))
fwrite(ob_timings, file="paper/ob_timings.csv")
summarised_timings_combined <- ob_timings[, .(chance50=quantile(ecdf(breach_times), probs=0.5), chance95=quantile(ecdf(breach_times), probs=0.95)), 
                                                    by=list(R0, VE, coverage)]

fwrite(summarised_timings_combined, file="paper/ob_timings_summary.csv")

ggplot(summarised_timings_combined, aes(x=coverage, y=chance50, colour=as.factor(R0))) + 
  geom_line() + 
  #geom_ribbon(aes(xmin=timelower95, xmax=timeupper95), colour=NA, alpha=0.1) +
  facet_wrap(~VE, scales="free")
