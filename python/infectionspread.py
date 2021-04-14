import numpy as np
import statistics   
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import nbinom

TRANSMISSION_POTENTIAL_MEAN = 1.3
TRANSMISSION_POTENTIAL_SD = 0.15
NUM_ESCAPED = 10

np.random.seed(1337)

def get_secondary_infections(transmission_potential, dispersion=0.1):

    # Python uses a messed up definition of a negative binomial, where you give the number of successes, P(success)
    # and get back the number of failures. The mean of this distribution is 1/p - 1 because there is 
    # E[trials to reach success] = 1/p, but the last trial is a success, so we subtract that one.
    # Converting that to transmission potential, if TP = 1.5, then the mean of the NB should be 1.5,
    # So 1.5 = 1/p - 1, i.e. p = 1/2.5.
    #secondary_infections = np.random.negative_binomial(1, 1.0/(transmission_potential+1))
    secondary_infections = nbinom.rvs(n=dispersion, p=dispersion/(dispersion+transmission_potential))
    # change p = 1/(1+tp) -> p = disp/(disp + tp)

    return (secondary_infections)

def get_infection_times(num_infections, scale=5.5, shape=1.8):

    #Numpy is a one parameter weibull so we just have to multiply:
    infection_times = scale*np.random.weibull(shape, num_infections)

    return (infection_times)

def run_generation(num_infected, infection_times, transmission_potential_mean = 1.3, transmission_potential_sd = 0.15, dispersion = 0.1, in_quarantine=False, in_isolation=False):
    secondary_infection_distribution = np.zeros(num_infected)
    next_generation_infections = 0
    next_generation_infections_time = []
    infection_days = 0
    for i in range(num_infected):
        #For reasons unclear to humans, the negative binomial, despite being an integer distribution, returns a float 
        #so we will have to cast this very shortly.
        transmission_potential = np.random.normal(transmission_potential_mean, transmission_potential_sd)
        while transmission_potential <= 0:
            transmission_potential = np.random.normal(transmission_potential_mean, transmission_potential_sd)    
        secondary_infection_distribution[i] = get_secondary_infections(transmission_potential, dispersion = dispersion)

        if in_quarantine:
            time_left_quarantine = np.random.poisson(1)
        else:
            time_left_quarantine = 0
        
        #Assuming they enter isolation 3 days after leaving quarantine.
        #Note if they leave quarantine on day 8, then isolate on day 11, the isolation is useless.
        if in_isolation:
            time_entered_isolation = 3+time_left_quarantine
        else:
            time_entered_isolation = np.inf

        if secondary_infection_distribution[i] > 0:
            next_infection_times = get_infection_times( int(secondary_infection_distribution[i]) )

            #Only keep infections that occurred after left quarantine
            #and before entered isolation
            infections_to_keep = (next_infection_times > time_left_quarantine) & (next_infection_times < time_entered_isolation)
            
            #Bit of casting magic to be able to iterate on booleans
            next_generation_infections += np.count_nonzero((infections_to_keep,))
            next_generation_infections_time += [ x+infection_times[i] for x in next_infection_times[infections_to_keep] ]
        
        #Assume infectious for 14 days, so they are in the community for 14-time spent in quarantine
        #Just make sure we don't falsely remove infection days
        #Assume 5.5 days is the infectious period.
        if time_left_quarantine <=5.5:
            #No benefit from isolation:
            if time_entered_isolation > 5.5:
                infection_days += 5.5 - time_left_quarantine
            else:
                #time_entered_isolation must be at least time_left_quarantine by above definition so this is always >0
                time_in_community = time_entered_isolation - time_left_quarantine
                infection_days += time_in_community
           
            

    return_dict = {
        "secondary_infection_distribution": secondary_infection_distribution,
        "infection_days": infection_days,
        "next_generation_infections": next_generation_infections,
        "infection_times": np.array(next_generation_infections_time).flatten(),
    }

    return (return_dict)
    

def run_outbreak(num_escaped, max_generations, generation_zero_times = None, max_infections = np.inf, transmission_potential_mean = 1.3, dispersion = 0.1, in_isolation = None, verbose=False):
    if generation_zero_times is None:
        generation_zero_times = np.zeros(num_escaped)
    
    if in_isolation is None:
        generation_zero = run_generation(num_escaped, generation_zero_times, dispersion = dispersion, in_isolation = False, in_quarantine=True, transmission_potential_mean=transmission_potential_mean)
    else:
        generation_zero = run_generation(num_escaped, generation_zero_times, dispersion = dispersion, in_isolation = True, in_quarantine=True, transmission_potential_mean=transmission_potential_mean)

    if verbose:
        print(num_escaped, "individuals escaped quarantine while positive.")
        print("Of those,", sum(generation_zero["secondary_infection_distribution"] > 0), "generated >0 infections (total of", int(sum(generation_zero["secondary_infection_distribution"])),"infections)")
        print("This totalled", generation_zero["infection_days"], "infectious days in the community and", generation_zero["next_generation_infections"], "infections seeded")

    iter = 0
    generations = []
    generations.append(generation_zero)
    while iter < max_generations and generations[iter]["next_generation_infections"] > 0 \
    and get_finalsize(generations) < max_infections:
        if in_isolation == "all":
            next_generation = run_generation(generations[iter]["next_generation_infections"], generations[iter]["infection_times"], dispersion = dispersion, in_isolation=True,
            transmission_potential_mean=transmission_potential_mean)
        else:
            next_generation = run_generation(generations[iter]["next_generation_infections"], generations[iter]["infection_times"], dispersion = dispersion, in_isolation=False,
            transmission_potential_mean=transmission_potential_mean)
        generations.append(next_generation)

        iter += 1
    
    return (generations)

def get_finalsize(outbreak):
    s = pd.Series( (x["next_generation_infections"] for x in outbreak))

    finalsize = sum(s)

    return (finalsize)

#Critical number of seedings:
def get_seedings(in_isolation, max_infected=50, max_seed_infected=30, iterations=10000, transmission_potential_mean=1.3, dispersion=0.1):
    died_out = np.zeros([max_seed_infected, max_infected+2])
    num_generations = np.zeros([max_seed_infected, 1000])
    categories = np.zeros([max_seed_infected, 3])
    infection_days = np.zeros([max_seed_infected, iterations])
    finalsize_ret = np.zeros([max_seed_infected, iterations])
    for initial_infected in tqdm(range(max_seed_infected)):
        
        for iter in range(iterations):
            #Range starts at 0:
            generation_zero_times = np.random.uniform(0, 12*7, size = initial_infected)
            outbreak = run_outbreak(initial_infected, max_infected+10, generation_zero_times, in_isolation=in_isolation, max_infections = max_infected, transmission_potential_mean = transmission_potential_mean, dispersion=dispersion)
            finalsize = get_finalsize(outbreak)
            if finalsize > max_infected:
                finalsize = max_infected
            died_out[int(initial_infected)][int(finalsize)] += 1
            num_generations[int(initial_infected)][len(outbreak)] += 1
            infection_days[int(initial_infected)][iter] = outbreak[0]["infection_days"]
            finalsize_ret[int(initial_infected)][iter] = finalsize

            #Category Green: No outbreak (i.e. 1 generation)
            if finalsize == 0:
                categories[initial_infected][0] += 1
            #Category yelllow: Finalsize < threshold
            elif finalsize < max_infected:
                categories[initial_infected][1] += 1
            #Categrory red: Finalsize >= threshold
            else:
                categories[initial_infected][2] += 1  

    return ({
        "died_out": died_out,
        "num_generations": num_generations,
        "categories": categories,
        "infection_days": infection_days,
        "finalsize": finalsize_ret
    })  
