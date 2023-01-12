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
m = 5
k_facs = 10

cover_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover") if run.endswith('.json') and "albe" in run]
supplier_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"supplier") if run.endswith('.json') and "albe" in run]
baseline_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"baseline") if run.endswith('.json') and "albe" in run]
flow_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"flow") if run.endswith('.json') and "albe" in run]

seeded_runs = {}

for file in cover_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover"/file) as f:
        result = json.load(f)
        seed = result[0]
        facs = result[1][str(k_facs)]["paths"]
        seeded_runs[seed] = {"cover": facs}

for file in supplier_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"supplier"/file) as f:
        result = json.load(f)
        seed = result[0]
        facs = result[1][str(k_facs)]["paths"]
        seeded_runs[seed]["supplier"] = facs

for file in baseline_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"baseline"/file) as f:
        result = json.load(f)
        seed = result[0]
        facs = result[1][str(k_facs)]["paths"]
        seeded_runs[seed]["baseline"] = facs

for file in flow_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"flow"/file) as f:
        result = json.load(f)
        seed = result[0]
        facs = result[1][str(k_facs)]["paths"]
        seeded_runs[seed]["flow"] = facs

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

percentile_list = [i/10 for i in range(800, 1001)]
#percentile_list = [i for i in range(1, 101)]

output_dict = {"cover": [], "supplier": [], "baseline":[], "flow":[], "percentile_list": percentile_list}
for seed, facs in seeded_runs.items():
    
    potential_facilities, location_directory, pid_assignment = data.get_data("albemarle", day, start, end, random_state=seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    radius_paths_cover = seeded_runs[seed]["cover"]
    radius_paths_supplier = seeded_runs[seed]["supplier"]
    radius_paths_baseline = seeded_runs[seed]["baseline"]
    radius_paths_flow = seeded_runs[seed]["flow"]
    
    objective_list_cover = calculate_objective_percentile_list(radius_paths_cover, location_directory, pid_assignment, percentile_list)
    objective_list_supplier = calculate_objective_percentile_list(radius_paths_supplier, location_directory, pid_assignment, percentile_list)
    objective_list_baseline = calculate_objective_percentile_list(radius_paths_baseline, location_directory, pid_assignment, percentile_list)
    objective_list_flow = calculate_objective_percentile_list(radius_paths_flow, location_directory, pid_assignment, percentile_list)
    
    output_dict["cover"].append([seed, {"paths": radius_paths_cover, "obj_list": objective_list_cover}])
    output_dict["supplier"].append([seed, {"paths": radius_paths_supplier, "obj_list": objective_list_supplier}])
    output_dict["baseline"].append([seed, {"paths": radius_paths_baseline, "obj_list": objective_list_baseline}])
    output_dict["flow"].append([seed, {"paths": radius_paths_flow, "obj_list": objective_list_flow}])
    

with open(PROJECT_ROOT/"Experiments"/"output"/f"albe_percentile_sensitivity_{m}_{k_facs}_narrow_robust.json", "w") as f:
    json.dump(output_dict, f)
