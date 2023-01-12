from FullyMobile import data
from FullyMobile.algorithm import *
from FullyMobile.utils import *
from FullyMobile import PROJECT_ROOT
import json
import os
import time
import random

day = 5
start = 6
end = 20

random_seeds = []

cover_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"cover") if run.endswith('.json') and "cville" in run]
prev_baseline_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"baseline") if run.endswith('.json') and "cville" in run]

for file in cover_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"cover"/file) as f:
        seed = json.load(f)["seed"]
        random_seeds.append(seed)

random_seeds = sorted(random_seeds)

for seed in random_seeds:
    
    if f"cville_movement_budget_long_{seed}.json" in prev_baseline_runs:
        print(seed)
        continue

    potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", day, start, end, random_state=seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    baseline_runs = {}
    
    #for travel in range(10, 55, 5):
    #for travel in range(55, 65, 5):
    #for travel in range(75, 85, 5):   
    
    for travel in range(10, 85, 5):
        
        m = travel/10
        
        baseline_runs[m] = {}
        
        for k_facs in range(3, 11):
            
            radius_paths = iterative_supplier(potential_facilities, location_directory, pid_assignment, k_facs, m, upper_limit = 15)
            #cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.01, upper_limit = 8, n_jobs = 40)
            objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
            print(m, k_facs, objective)
            baseline_runs[m][k_facs] = {"paths": radius_paths, "objective": objective}
    
    with open(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"baseline"/f"cville_movement_budget_long_{seed}.json", "w") as f:
        json.dump({"seed": seed, "baseline_runs": baseline_runs}, f)