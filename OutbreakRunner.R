library(reticulate)
library(data.table)
library(dplyr)
library(ggplot2)
library(purrr)

#Comment this out if you are not using a custom python environment
use_condaenv("./python/envs")

source_python("python/infectionspread.py")

#' A thin wrapper around the python function run_outbreak contained in infectionspread.py
#' Default stopping condition is when there are 0 new cases *which may never be met*.
#' 
#' Arguments:
#' initial_infected: Number of initial infections to seed
#' transmission_potential_mean: Mean transmission potential
#' dispersion: Dispersion. Default 0.1
#' max_infections: Maximum number of infections, beyond which the simulation will terminate. Default Inf
#' max_generations: Maximum number of branching process generations, beyond which the simulation will terminate. Default 25. *Not currently affecting anything*
#' max_time: Maximum simulation time
run_outbreak.wrapper <- function(initial_infected, transmission_potential_mean, 
                                 dispersion = 0.1, max_infections = Inf, max_generations = 25, max_time = 100) {
  
  run_outbreak(as.integer(initial_infected), as.integer(max_generations), 
               transmission_potential_mean = transmission_potential_mean, 
               dispersion = dispersion,
               max_infections = max_infections,
               max_time = max_time)
}

#' Counts the number of new infections and cumulative number of infections, by day, in a given outbreak.
#' Arguments:
#' outbreak: Output from run_outbreak.wrapper
#' max_time: Maximum time, beyond which results are trunacted. Required due to simulation artifact where "long" chains will appear decreasing
#' but are just because other chains have stopped.
get_infections <- function(outbreak, max_time) {
  infection_times <- map(outbreak, "infection_times") %>% do.call(c, .)
  if (is.null(infection_times)) {
    infections <- data.table(day = 0, cumulative_infected = 0, new_infections = NA)
    
    return (infections)
  }
  
  infections <- data.table(day = 0:max_time, cumulative_infected = sapply(0:max_time, function(i) { sum(infection_times < i) }) )
  infections[, new_infections := c(NA, diff(cumulative_infected))]

  return (infections)  
}


#' Wraps up run_outbreak.wrapper and get_infections to allow for 
#' convenient mass simulation of outbreaks.
#' #' Default stopping condition is when there are 0 new cases *which may never be met*.
#' 
#' Arguments:
#' initial_infected: Number of initial infections to seed
#' transmission_potential_mean: Mean transmission potential
#' iterations: Number of replications to perform
#' dispersion: Dispersion. Default 0.1
#' max_infections: Maximum number of infections, beyond which the simulation will terminate. Default Inf
#' max_generations: Maximum number of branching process generations, beyond which the simulation will terminate. Default 25.
run_outbreak_multiple <- function(initial_infected, transmission_potential_mean, iterations, 
                                  dispersion = 0.1, max_infections = Inf, max_generations = 25, max_time = 100) {
  
  lapply(1:iterations, function(iter) {
    
    ob <- run_outbreak.wrapper(initial_infected = initial_infected, 
                               transmission_potential_mean = transmission_potential_mean, 
                               dispersion = dispersion, 
                               max_infections = max_infections,
                               max_generations = max_generations,
                               max_time = max_time)
    
    infections <- get_infections(ob, max_time)
    infections[, iteration := iter]
  }) %>% rbindlist()
  
}
