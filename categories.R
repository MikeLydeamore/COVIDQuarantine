library(reticulate)
library(data.table)
library(dplyr)
library(ggplot2)

use_condaenv("./python/envs")

source_python("python/infectionspread.py")
py_run_string("from tqdm.notebook import tqdm
from matplotlib import gridspec")

#Uncomment seedings_isolated for the model where individuals are isolated 3 days into the GI (cutting off all infection times > 3)
#This block often runs strangely so you'll need to highlight the whole thing and run it.
py_run_string("transmission_potentials = [1.1, 1.5, 2]

iterations=10000
seedings = []
seedings_isolated = []

for TRANSMISSION_POTENTIAL_MEAN in tqdm(transmission_potentials):
    seedings.append(get_seedings(in_isolation=None, max_infected=20, max_seed_infected=30, transmission_potential_mean=TRANSMISSION_POTENTIAL_MEAN, iterations=iterations, dispersion = 1.0))
    #seedings_isolated.append(get_seedings(in_isolation='all', max_infected=20, max_seed_infected=30, transmission_potential_mean=TRANSMISSION_POTENTIAL_MEAN, iterations=iterations))"
)


categories <- lapply(1:3, function(i) {
  categories <- as.data.table(py$seedings[[i]]$categories)
  
  colnames(categories) <- c("zero","middle","top")
  categories[, initial_infected := 0:(.N-1)][, transmission_potential := unlist(py$transmission_potentials)[i]]
  
  cat_noisolation <- melt(categories, id.vars=c("initial_infected","transmission_potential"), variable.name = "category")
  cat_noisolation[, isolation := "none"]
  
  #categories <- as.data.table(py$seedings_isolated[[i]]$categories)
  
  #colnames(categories) <- c("zero","middle","top")
  #categories[, initial_infected := 0:(.N-1)][, transmission_potential := unlist(py$transmission_potentials)[i]]
  
  #cat_isolation <- melt(categories, id.vars=c("initial_infected","transmission_potential"), variable.name = "category")
  #cat_isolation[, isolation := "all"]
  
  return (cat_noisolation)
  
}) %>% rbindlist()

categories$category <- factor(categories$category, levels=c("zero","middle","top"), labels = c("Zero cases", "<20 cases", ">20 cases"))
categories$isolation <- factor(categories$isolation, levels=c("none","all"), label=c("No isolation", "Isolation after 3 days"))
ggplot(categories[initial_infected != 0 & isolation == "No isolation",], aes(x=initial_infected*4.5, y=value/10000, fill=category)) + geom_col(position = position_stack(reverse=F), width=4.5) +
  facet_wrap(~transmission_potential, ncol=1) +
  labs(x="First generation person-days infectious in the community", y="Probality of category") +
  theme_bw() + 
  theme(legend.title = element_blank(), panel.spacing.y = unit(1, "lines")) +
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("#B3CDE3", "#8C96C6", "#810F7C"))



probability_arrivals <- function(arrivals, categories) {
  #2.5 infection-days per 10,000 arrivals.
  infection_days_per_week <- arrivals/10000 * 2.5
  #Weeks to get 4.5 infection-days:
  weeks_to_one_case <- 4.5/infection_days_per_week

  #Number of people in community to have 95% chance of large outbreak:
  people_required <- categories[isolation == "No isolation" & category == ">20 cases" & value > 0.95*10000, .(people_required=min(initial_infected)), by=transmission_potential]
  
  #Don't want to edit the base frame:
  df <- categories[isolation == "No isolation" & category == ">20 cases",]
  
  people_required[, weeks_to_outbreak := people_required*weeks_to_one_case]
  
  df[, prob := value/10000][, time := today()+weeks(initial_infected*weeks_to_one_case)]
  
}


plot_frame <- lapply(c(3000), function(arrivals) {
  
  x <- probability_arrivals(arrivals = arrivals, categories)
  x$cap <- arrivals
  
  x1 <- probability_arrivals(arrivals = 0.5*arrivals, categories)
  x1$cap <- 0.5*arrivals
  
  x2 <- probability_arrivals(arrivals = 2*arrivals, categories)
  x2$cap <- 2*arrivals
  
  ret <- rbindlist(list(x,x1,x2))
  
  ret$arrivals <- arrivals
  
  return (ret)
}) %>% rbindlist()


facet_lookup <- c("1.1"="Transmission potential: 1.1",
                  "1.5"="Transmission potential: 1.5",
                  "2"="Transmission potential: 2.0")

ggplot(plot_frame, 
       aes(y=prob)) + 
  geom_line(aes(x=time, colour=as.factor(cap)), size=1) + 
  #geom_area(aes(x=today()+weeks(initial_infected*weeks_to_one_case)), fill="#810F7C", alpha=0.5) +
  #geom_line(aes(x=today()+weeks(2*initial_infected*weeks_to_one_case)), linetype="dashed", size=1) + 
  #geom_area(aes(x=today()+weeks(2*initial_infected*weeks_to_one_case)), fill="#B3CDE3", alpha=0.5) +
  facet_wrap(~transmission_potential, ncol=1, scales="free", labeller = as_labeller(facet_lookup)) +
  labs(x="Date", y="Probability of major outbreak") +
  cowplot::theme_cowplot() + 
  labs(colour="Arrival cap (per week)") +
  theme(panel.spacing.y = unit(1, "lines")) +
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_date(expand=c(0,0), breaks = ymd("2021-03-01")+months(0:11), date_labels="%e/%m") +
  scale_fill_manual() +
  coord_cartesian(xlim=(c(today(), ymd("2021-12-31")))) +
  theme(strip.background = element_rect(fill = NA, colour = "black"), legend.position = "bottom", legend.justification = "center")

ggsave("outbreak_likelihoods_dispersion08.png", plot = last_plot(), type="cairo", width = 7.5, height = 10, units="in")
