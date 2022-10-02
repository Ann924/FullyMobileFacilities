from FullyMobile import data
from FullyMobile.algorithm import *
from FullyMobile.utils import *
from FullyMobile import PROJECT_ROOT
import json
import os
import time

day = 5
start = 6
end = 20

potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", day, start, end)
print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))

cover_runs = {}

for travel in [1, 3, 5]:
    
    m = travel
    
    cover_runs[m] = {}
    
    for k_facs in [3, 5, 10]:
        opt_radius, radius_paths = cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.01, upper_limit = 8, n_jobs = 40)
        objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
        print(m, objective)
        cover_runs[m][k_facs] = {"paths": radius_paths, "objective": objective}

with open(PROJECT_ROOT/"Experiments"/"output"/"cville_movement_budget.json", "w") as f:
    json.dump(cover_runs, f)