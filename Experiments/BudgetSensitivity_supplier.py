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
for k_facs in range(3, 11):
    radius_paths = stationary_supplier(potential_facilities, location_directory, pid_assignment, k_facs)
    objective = calculate_assignment_cost(radius_paths, location_directory, pid_assignment)
    print(objective)
    cover_runs[k_facs] = {"paths": radius_paths, "objective": objective}

with open(PROJECT_ROOT/"Experiments"/"output"/"cville_budget_sensitivity_supplier.json", "w") as f:
    json.dump(cover_runs, f)