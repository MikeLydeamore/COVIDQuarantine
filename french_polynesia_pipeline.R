library(reticulate)
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)

use_condaenv("./python/envs")
source_python("./python/quarantine_pipeline.py")
source_python("./python/infectionspread.py")
source("quarantine_pipeline_functions.R")


scenarios <- list.files("D:/The University of Melbourne/Cameron Zachreson - PYF_Scenarios_for_Michael/", pattern = "Run*")

tps <- fread("data/penetration/tps_ML_model.csv")

combinations <- expand.grid(scenarios, tps$new_tp_delta)
colnames(combinations) <- c("scenario", "tp")

ob_probs <- mapply(SIMPLIFY = FALSE, scenario=combinations$scenario, tp=combinations$tp, FUN=function(scenario, tp) {
  message(paste0("Running scenario ", scenario, " at TP=", tp))
  base_fname <- paste0("D:/The University of Melbourne/Cameron Zachreson - PYF_Scenarios_for_Michael/", scenario)
  
  gen0tp <- tps[new_tp_delta == tp, new_tp_import_unvac_delta]
  gen0tp_vaccinated <- tps[new_tp_delta == tp, ifelse(scenario %like% "AZ2|COM2", new_tp_import_vac_pf_delta, new_tp_import_vac_az1_delta)]
  
  ob_prob_traveller <- probability_of_major_outbreak(100000, tp, quarantine_type = "HomeQuar", type="traveller",
  travellers_vaccinated=NA, filename=paste0(base_fname, "/Traveller_breach.csv"),
  tpgen0=gen0tp, tpgen0_vaccinated=gen0tp_vaccinated)
  #ob_prob_worker <- probability_of_major_outbreak(1000, 1.7, quarantine_type = "HomeQuar", type="worker",
  #                                                  travellers_vaccinated=NA, filename=paste0(base_fname, "/Worker_breach.csv"))
  ret <- data.table(scenario, ob_prob_traveller, tp)
  return (ret)
})
ob_probs <- rbindlist(ob_probs)


ob_timings <- mapply(SIMPLIFY=FALSE, 
                     scenario=combinations$scenario, tp=combinations$tp,
                     FUN=function(scenario, tp) {
                       message(paste0("Running scenario ", scenario, " at TP=", tp))
                       base_fname <- paste0("D:/The University of Melbourne/Cameron Zachreson - PYF_Scenarios_for_Michael/", scenario)
                       ob_prob <- ob_probs[scenario == scenario & tp == tp, ob_prob_traveller]
                       summary_file <- readLines(paste0(base_fname, "/",substr(scenario, 1, nchar(as.character(scenario))-11), "__results_summary.txt"))
                       travellers_discharged_line <- summary_file[which(summary_file %like% "travellers discharged per week")[1]]
                       
                       travellers_per_week <- as.numeric(regmatches(travellers_discharged_line, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",travellers_discharged_line)))
                       scaling_factor <- round(2000 / travellers_per_week)
                       
                       
                       breach_times_outbreak <- timing_of_major_outbreak(1e6, ob_prob=ob_prob, 
                                                                         quarantine_type="HomeQuar", quarantine_duration=10, 
                                                                         travellers_vaccinated="Vacc",
                                                                         infected_proportion=0.01,
                                                                         num_arrivals=2000,
                                                                         scaling_factor=scaling_factor,
                                                                         filename=paste0(base_fname, "/Traveller_breach.csv"))
                       
                       data.table(scenario, tp, ob_timings = breach_times_outbreak)
                     }) %>% rbindlist()


ob_timings[, .(chance50=quantile(ecdf(ob_timings), probs=0.5), chance95=quantile(ecdf(ob_timings), probs=0.95)), 
  by=list(scenario, tp)][order(scenario, tp)]

saveRDS(ob_timings, "D:/The University of Melbourne/Cameron Zachreson - PYF_Scenarios_for_Michael/outbreak_timings.rds")
fwrite(ob_timings[, .(chance50=quantile(ecdf(ob_timings), probs=0.5), chance95=quantile(ecdf(ob_timings), probs=0.95)), 
                  by=list(scenario, tp)][order(scenario, tp)], file="D:/The University of Melbourne/Cameron Zachreson - PYF_Scenarios_for_Michael/outbreak_timings_summary.csv")

fwrite(ob_probs, "D:/The University of Melbourne/Cameron Zachreson - PYF_Scenarios_for_Michael/outbreak_probabilities.csv")
