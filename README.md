# Branching process model to simulate outbreaks in the community

This repository contains the work to take "breach events" (generated by a separate model) and seed them in the community which has a defined level of vaccine coverage (with variable efficacy).

The file [`paper_pipeline.R`](paper_pipeline.R) outlines the process to run the model. The statistics are generated in two steps:

1. Calculate the probability of any given breach event resulting in an outbreak (defined as >5 cases that occur before eradication).
2. Determine the timing of outbreaks, based on the frequency distribution of breach events.

These two steps allow for calculating the probability of an outbreak occuring in some number of days, or equivalently the expected number of outbreaks over the year.

# Model code
The underlying model is written in python (see in directory [`python`](python/), and is imported into R using `reticulate`. 

Breaches resulting in local transmission are generally rare events, so it is recommended to use a very large number of iterations for step 1 of the above algorithm, at least 100,000.
