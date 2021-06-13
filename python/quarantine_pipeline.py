import pandas as pd
import random
import numpy as np

sys.path.append("./python")

import infectionspread

def calculate_concurrent_outbreaks(transmission_potential, breaches_over_horizon, quarantine_type, quarantine_duration, df_traveller):
    
    time_between_breaches = df_traveller["time_discharged"][1:].values - df_traveller["time_discharged"][0:-1].values
    breach_times = np.cumsum(random.choices(time_between_breaches, k=breaches_over_horizon))
    gen_zero_times = random.choices(df_traveller["time_left_quarantine"], k=breaches_over_horizon)
    ob = []

    for iter in range(breaches_over_horizon):
    
        ob.append(run_outbreak(num_escaped=1, max_generations=np.inf, max_infections = 20,
                               generation_zero_time_left_quarantine=(gen_zero_times[iter],),
                               transmission_potential_mean=transmission_potential, max_time=50 ) )

    results = {
    'finalsize': [get_finalsize(outbreak) for outbreak in ob],
    'breach_times': breach_times,
    'days_in_community': gen_zero_times
    }

    df = pd.DataFrame(data=results)

    df["outbreak_end_time"]=df["breach_times"]+21
    df["is_outbreak"]=df["finalsize"]>=5

    df_outbreaks = df[df.is_outbreak].sort_values(by="breach_times")
    
    max_concurrent_outbreaks = 0
    if len(df_outbreaks) > 0:
        for i in range(len(df_outbreaks)-1):
            concurrent_outbreaks = 1
            end_time = df_outbreaks["outbreak_end_time"].values[i]
            for j in range(i+1, len(df_outbreaks)):
                start_time = df_outbreaks["breach_times"].values[j]
                if start_time < end_time:
                    concurrent_outbreaks += 1
            max_concurrent_outbreaks = max(max_concurrent_outbreaks, concurrent_outbreaks)

    return ({
        'max_concurrent_outbreaks': max_concurrent_outbreaks,
        'df': df
    })
    

def concurrent_outbreak_probabilities(transmission_potential, breaches_over_horizon, quarantine_type, quarantine_duration, max_iterations=10000):
    #Data loading  block:
    filename = "..\data\penetration\TUnvacc_%s_%sd_traveller_breach_timeseries.csv" % (quarantine_type, str(quarantine_duration))
    df_traveller = pd.read_csv(filename)
    df_traveller["time_left_quarantine"] = df_traveller["t_incubation"] + df_traveller["t_post_incubation"] - df_traveller["days_infectious_community"]

    if breaches_over_horizon is None:
        num_breaches = len(df_traveller)
        time_period = max(df_traveller["time_discharged"])
        avg_breaches = num_breaches / time_period

        time_horizon = 365
        breaches_over_horizon = round(avg_breaches * time_horizon)

        print(str(breaches_over_horizon)+" breaches for TP="+str(transmission_potential)+" and quarantine length="+str(quarantine_duration))

    concurrent_outbreaks = []
    for i in range(max_iterations):
        x=calculate_concurrent_outbreaks(transmission_potential=transmission_potential, breaches_over_horizon=breaches_over_horizon, 
                                         quarantine_type=quarantine_type, quarantine_duration=quarantine_duration, df_traveller=df_traveller)
        concurrent_outbreaks.append(x["max_concurrent_outbreaks"])

    ret = pd.Series(concurrent_outbreaks).value_counts(sort=False) / max_iterations

    return (ret)
