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

m = 5

cover_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover") if run.endswith('.json') and "albe" in run]

prev_flow_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"flow") if run.endswith('.json') and "albe" in run]

random_seeds = []

for file in cover_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover"/file) as f:
        seed = json.load(f)[0]
        random_seeds.append(seed)

for seed in random_seeds:
    
    if f"albe_budget_sensitivity_flow_{seed}.json" in prev_flow_runs:
        print(seed)
        continue
    
    potential_facilities, location_directory, pid_assignment = data.get_data("albemarle", day, start, end, random_state=seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    flow_runs = {}
    for k_facs in range(6, 21):
        radius_paths = flow_search(potential_facilities, location_directory, pid_assignment, k_facs, m, upper_limit = 20)
        objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
        print(objective)
        flow_runs[k_facs] = {"paths": radius_paths, "objective": objective}
    
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"flow"/f"albe_budget_sensitivity_flow_{seed}.json", "w") as f:
        json.dump([seed, flow_runs], f)