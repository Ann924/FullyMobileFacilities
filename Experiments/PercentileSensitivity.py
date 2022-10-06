from FullyMobile import data
from FullyMobile.algorithm import *
from FullyMobile.utils import *
from FullyMobile import PROJECT_ROOT
import json
import math
import os
import time

day = 5
start = 6
end = 20
m = 3
k_facs = 5

potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", day, start, end)
print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))

def calculate_objective_percentile_list(facilities, location_directory, pid_assignment, percentile_list):
    distance_list = list()
    
    START = min(pid_assignment.keys())
    END = max(pid_assignment.keys())
    
    for hr in range(START, END):
        for pid, client in pid_assignment[hr]:
            coord_client = (location_directory[client]['latitude'], location_directory[client]['longitude'])
            client_cost = min([calculate_dist(coord_client, (location_directory[fac[hr-START]]['latitude'], location_directory[fac[hr-START]]['longitude'])) for fac in facilities])
            distance_list.append(client_cost)
    
    distance_list = sorted(distance_list)
    return [distance_list[math.ceil(len(distance_list)*p/100)-1] for p in percentile_list if p!=0]

percentile_list = [i for i in range(0, 101)]

opt_radius_cover, radius_paths_cover = cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.1, upper_limit = 4, n_jobs = 40)
radius_paths_supplier = stationary_supplier(potential_facilities, location_directory, pid_assignment, k_facs)

objective_list_cover = calculate_objective_percentile_list(radius_paths_cover, location_directory, pid_assignment, percentile_list)
objective_list_supplier = calculate_objective_percentile_list(radius_paths_supplier, location_directory, pid_assignment, percentile_list)

with open(PROJECT_ROOT/"Experiments"/"output"/f"percentile_sensitivity_{m}_{k_facs}.json", "w") as f:
    data = {
        "cover":{
            "paths": radius_paths_cover,
            "obj_list": objective_list_cover
        },
        "supplier":{
            "paths": radius_paths_cover,
            "obj_list": objective_list_supplier
        },
        "percentile_list": percentile_list
    }
    json.dump(data, f)