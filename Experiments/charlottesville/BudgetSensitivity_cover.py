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

for i in range(0, 4):
    seed = random.randint(0, 1e9)
    
    potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", day, start, end, random_state=seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    m = 5
    
    cover_runs = {}
    for k_facs in range(3, 11):
        opt_radius, radius_paths = cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.1, upper_limit = 4, n_jobs = 40)
        objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
        print(objective)
        cover_runs[k_facs] = {"paths": radius_paths, "objective": objective}
    
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover"/f"cville_budget_sensitivity_cover_{seed}.json", "w") as f:
        json.dump([seed, cover_runs], f)