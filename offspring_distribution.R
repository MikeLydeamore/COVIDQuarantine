library(stats)
library(dplyr)
library(data.table)

offspring <- lapply(c(1.1, 1.5, 2), function(tp) {
  data.table(offspring=rnbinom(n = 10000, size=0.1, mu = rnorm(n=10000, mean=tp, sd=0.15)), transmission_potential=tp)
}) %>% rbindlist()

facet_lookup <- c("1.1"="Transmission potential: 1.1",
                  "1.5"="Transmission potential: 1.5",
                  "2"="Transmission potential: 2.0")

ggplot(offspring, aes(x=offspring, y=after_stat(density))) + 
  facet_wrap(~transmission_potential, labeller = as_labeller(facet_lookup), ncol=1, scales="free_x") + 
  geom_histogram(binwidth=1, fill='#B3CDE3', colour="black") +
  labs(x="Number of secondary infections", y="Probability") +
  cowplot::theme_cowplot() +
  theme(strip.background = element_rect(fill = NA, colour = "black"), legend.position = "bottom", legend.justification = "center") +
  scale_x_continuous(limits=c(-1, 20), expand=c(0,0))
  
ggsave("offspring_distribution.png", plot = last_plot(), type="cairo", width = 7.5, height = 10, units="in")
