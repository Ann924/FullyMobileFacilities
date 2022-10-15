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

cover_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover") if run.endswith('.json') and "albe" in run]
random_seeds = []

for file in cover_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover"/file) as f:
        seed = json.load(f)[0]
        random_seeds.append(seed)

for seed in random_seeds:
    
    potential_facilities, location_directory, pid_assignment = data.get_data("albemarle", day, start, end, random_state=seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    cover_runs = {}
    for k_facs in range(6, 21):
        radius_paths = stationary_supplier(potential_facilities, location_directory, pid_assignment, k_facs, lower_limit = 0.1, upper_limit = 20)
        objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
        print(objective)
        cover_runs[k_facs] = {"paths": radius_paths, "objective": objective}
    
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"supplier"/f"albe_budget_sensitivity_supplier_{seed}.json", "w") as f:
        json.dump([seed, cover_runs], f)