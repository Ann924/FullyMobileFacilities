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

for file in cover_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"cover"/file) as f:
        seed = json.load(f)["seed"]
        random_seeds.append(seed)

random_seeds = sorted(random_seeds)

for seed in random_seeds[::-1]:
    
    prev_flow_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"flow") if run.endswith('.json') and "cville" in run]
    
    if f"cville_movement_budget_long_{seed}.json" in prev_flow_runs:
        print(seed)
        continue

    potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", day, start, end, random_state=seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    flow_runs = {}
    
    #for travel in range(10, 55, 5):
    #for travel in range(55, 65, 5):
    #for travel in range(75, 85, 5):
    
    seed_completed = False
    
    for travel in range(10, 85, 5):
        
        m = travel/10
        
        flow_runs[m] = {}
        
        for k_facs in range(3, 11):
            
            prev_flow_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"flow") if run.endswith('.json') and "cville" in run]
    
            if f"cville_movement_budget_long_{seed}.json" in prev_flow_runs:
                seed_completed = True
                print(seed)
                break
            
            radius_paths = flow_search(potential_facilities, location_directory, pid_assignment, k_facs, m, upper_limit=15)
            #radius_paths = iterative_supplier(potential_facilities, location_directory, pid_assignment, k_facs, m, upper_limit = 15)
            #cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.01, upper_limit = 8, n_jobs = 40)
            objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
            print(seed, m, k_facs, objective)
            flow_runs[m][k_facs] = {"paths": radius_paths, "objective": objective}
        
        if seed_completed:
            break
    
    if seed_completed:
        continue
    
    with open(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/"flow"/f"cville_movement_budget_long_{seed}.json", "w") as f:
        json.dump({"seed": seed, "flow_runs": flow_runs}, f)