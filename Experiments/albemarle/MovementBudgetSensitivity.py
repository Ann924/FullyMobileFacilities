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

seed = random.randint(0, 1e9)

potential_facilities, location_directory, pid_assignment = data.get_data("albemarle", day, start, end, random_state=seed)
print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))

cover_runs = {}

#for travel in range(75, 105, 5):
#for travel in range(55, 75, 5):
#for travel in range(10, 55, 5):

for travel in range(10, 105, 5):
    
    m = travel/10
    
    cover_runs[m] = {}
    
    #[i for i in range(8, 20, 2) if i not in [6, 10, 20]]
    for k_facs in range(6, 21, 2):
        opt_radius, radius_paths = cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.01, upper_limit = 20, n_jobs = 40)
        objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
        print(m, k_facs, objective)
        cover_runs[m][k_facs] = {"paths": radius_paths, "objective": objective}

    # keep inside movement for-loop so partial plots can still be created
    with open(PROJECT_ROOT/"Experiments"/"output"/"movement_sensitivity"/f"albe_movement_budget_long_{seed}.json", "w") as f:
        json.dump({"seed": seed, "cover_runs": cover_runs}, f)