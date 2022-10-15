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

range_dict = {
    6: (8, 25),
    7: (8, 25),
    8: (5, 20),
    9: (5, 20),
    10: (5, 15),
    11: (3, 10),
    12: (3, 10),
    13: (3, 10),
    14: (3, 10),
    15: (3, 10),
    16: (0.1, 8),
    17: (0.1, 8),
    18: (0.1, 8),
    19: (0.1, 8),
    20: (0.1, 8)
}

for i in range(0, 5):
    seed = random.randint(0, 1e9)
    #seed=42
    
    potential_facilities, location_directory, pid_assignment = data.get_data("albemarle", day, start, end, random_state=seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    m = 5
    
    cover_runs = {}
    for k_facs in range(6, 21):
        
        if k_facs in range_dict.keys():
            lower_limit = range_dict[k_facs][0]
            upper_limit = range_dict[k_facs][1]
        else:
            lower_limit = 0.1
            upper_limit = 25
        
        opt_radius, radius_paths = cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = lower_limit, upper_limit = upper_limit, n_jobs = 40)

        next_k = k_facs+1
        if next_k in range_dict.keys():
            range_dict[next_k] = (range_dict[next_k][0], opt_radius+0.5)

        objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
        print(i, k_facs, objective)
        cover_runs[k_facs] = {"paths": radius_paths, "objective": objective}
    
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover"/f"albe_budget_sensitivity_cover_{seed}.json", "w") as f:
        json.dump([seed, cover_runs], f)