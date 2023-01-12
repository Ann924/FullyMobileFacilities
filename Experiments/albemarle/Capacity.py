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

def calculate_capacity_redistribution(facilities, location_directory, pid_assignment, flexibility=0, constraint = 15):
    max_dist = 0
    
    START = min(pid_assignment.keys())
    END = max(pid_assignment.keys())
    
    for hr in range(START, END):
        clients = set(client[1] for client in pid_assignment[hr])
        print(len(clients))
        radius = redistribute([facilities[i][hr-START] for i in range(len(facilities))], location_directory, clients, flexibility=flexibility, constraint=constraint)
        
        if radius > max_dist:
            max_dist = radius
        print(hr, radius)
    
    return max_dist

def redistribute(facilities, location_directory, clients, flexibility=0, constraint = 15):
    facility_to_clients = {fac:[] for fac in facilities}
    client_facility_distances = {}
    
    max_dist = 0
    
    for c in clients:
        coord_client = (location_directory[c]['latitude'], location_directory[c]['longitude'])
        
        client_facility_distances[c] = {}
        
        min_fac = (10e9, -1)
        for fac in facilities:
            coord_fac = (location_directory[fac]['latitude'], location_directory[fac]['longitude'])
            distance = calculate_dist(coord_client, coord_fac)
            if distance < min_fac[0]:
                min_fac = (distance, fac)
            
            client_facility_distances[c][fac] = distance
        
        facility_to_clients[min_fac[1]].append(c)
        
        if min_fac[0] > max_dist:
            max_dist = min_fac[0]
    
    #constraint = math.ceil(len(clients)/len(facilities)) + flexibility
    #constraint = 15
    CAPACITY = {fac: constraint*facilities.count(fac) for fac in facilities}
    
    facility_counts = {fac: len(c) for fac, c in facility_to_clients.items()}
    facility_excess = {fac: facility_counts[fac]-CAPACITY[fac] for fac in facility_counts.keys()}
    
    #print(facility_counts, CAPACITY)
    
    #print(facility_counts, facility_excess, CAPACITY)
    
    largest_fac = max(facility_excess, key=facility_excess.get)
    
    counter = 0
    
    while facility_counts[largest_fac] > CAPACITY[largest_fac]:
        
        redist_choice = []
        for c in facility_to_clients[largest_fac]:
            options = [(val, key) for key, val in client_facility_distances[c].items() if key!=largest_fac and facility_excess[largest_fac]-facility_excess[key]>1]
            if len(options)>0:
                min_option = min(options)
                redist_choice.append((min_option[0], min_option[1], c))

        dist, fac, client = min(redist_choice)

        facility_to_clients[largest_fac].remove(client)
        facility_counts[largest_fac]-=1
        facility_excess[largest_fac]-=1

        facility_to_clients[fac].append(client)
        facility_counts[fac]+=1
        facility_excess[fac]+=1

        counter += 1
        
        largest_fac = max(facility_excess, key=facility_excess.get)
        
    print(facility_counts, CAPACITY)
    
    return max(client_facility_distances[c][fac] for fac, assigned_clients in facility_to_clients.items() for c in assigned_clients)


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


results = {}
for constraint in range(14, 31):
    
    output_dict = {"cover": [], "supplier": [], "baseline":[], "flow":[]}
    
    for seed, facs in seeded_runs.items():
        
        potential_facilities, location_directory, pid_assignment = data.get_data("albemarle", day, start, end, random_state=seed)
        print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
            
        radius_paths_cover = seeded_runs[seed]["cover"]
        radius_paths_supplier = seeded_runs[seed]["supplier"]
        radius_paths_baseline = seeded_runs[seed]["baseline"]
        radius_paths_flow = seeded_runs[seed]["flow"]
        
        objective_cover = calculate_capacity_redistribution(radius_paths_cover, location_directory, pid_assignment, constraint=constraint)
        objective_supplier = calculate_capacity_redistribution(radius_paths_supplier, location_directory, pid_assignment, constraint=constraint)
        objective_baseline = calculate_capacity_redistribution(radius_paths_baseline, location_directory, pid_assignment, constraint=constraint)
        objective_flow = calculate_capacity_redistribution(radius_paths_flow, location_directory, pid_assignment, constraint=constraint)
        
        output_dict["cover"].append([seed, {"paths": radius_paths_cover, "objective": calculate_assignment_cost(radius_paths_cover, location_directory, pid_assignment), "capacity_cost": objective_cover}])
        output_dict["supplier"].append([seed, {"paths": radius_paths_supplier, "objective": calculate_assignment_cost(radius_paths_supplier, location_directory, pid_assignment), "capacity_cost": objective_supplier}])
        output_dict["baseline"].append([seed, {"paths": radius_paths_baseline, "objective": calculate_assignment_cost(radius_paths_baseline, location_directory, pid_assignment), "capacity_cost": objective_baseline}])
        output_dict["flow"].append([seed, {"paths": radius_paths_flow, "objective": calculate_assignment_cost(radius_paths_flow, location_directory, pid_assignment), "capacity_cost": objective_flow}])
        
    results[constraint] = output_dict

with open(PROJECT_ROOT/"Experiments"/"output"/f"albe_capacity_{m}_{k_facs}_constraint.json", "w") as f:
    json.dump(results, f)
