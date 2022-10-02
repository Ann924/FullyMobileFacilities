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
k_facs = 10

cover_runs = {}

potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", day, start, end)
print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))

for travel in range(10, 80, 5):
    
    m = travel/10
    
    if os.path.exists(PROJECT_ROOT/f"reachable_{m}.json"):
        with open(PROJECT_ROOT/f"reachable_{m}.json") as f:
            reachable_dict = {int(k): v for k, v in json.load(f).items()}
    else:
        reachable_dict = calculate_reachable(potential_facilities, location_directory, 5)
        with open(PROJECT_ROOT/f"reachable_{m}.json", "w") as f:
            json.dump(reachable_dict, f)

    opt_radius, radius_paths = cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.01, upper_limit = 8, n_jobs = 40)
    objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
    print(m, objective)
    cover_runs[m] = {"paths": radius_paths, "objective": objective}
    
with open(PROJECT_ROOT/"Experiments"/"output"/"cville_movement_sensitivity.json", "w") as f:
    json.dump(cover_runs, f)