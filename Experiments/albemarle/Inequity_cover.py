from FullyMobile.data import *
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

def income_bracket_mapping(x, income_to_weight):
    
    for income, weight in income_to_weight.items():
        if x>=income[0] and x<income[1]:
            return weight

def generate_weight_mapping(county_name, income_to_weight):
    filename = PROJECT_ROOT/"Data"/f"usa_va_{county_name}_adult_activity_location_assignment_week.csv"
    df_county = pd.read_csv(filename)[['pid', 'hid']]

    filename = PROJECT_ROOT/"Data"/"va_household.csv"
    df = pd.read_csv(filename)
    df['income'] = df['hh_income']/df['hh_size']
    
    df['income_group'] = df['income'].apply(lambda x: income_bracket_mapping(x, income_to_weight))
    
    df_hid = df.set_index('hid')
    df_county = df_county.join(df_hid, on='hid')
    df_county = df_county[['hid', 'pid', 'income', 'income_group']].groupby(['pid']).mean()
    
    print([sum(df_county['income_group']==val) for val in income_to_weight.values()])
    
    return df_county['income_group'].to_dict()

def cover_pid_group(fac_loc, loc_assignment, location_directory, pid_weight_map, radius):
    
    cover_count = 0
    fac_coord = (location_directory[fac_loc]["latitude"], location_directory[fac_loc]["longitude"])
    
    coverage = set()
    
    for ind in range(len(loc_assignment)):
        loc = loc_assignment[ind]
        pid = loc[0]
        pid_weight = pid_weight_map[pid]
        lid = loc[1]
        
        lid_coord = (location_directory[lid]["latitude"], location_directory[lid]["longitude"])
        if calculate_dist(fac_coord, lid_coord) <= pid_weight * radius:
            cover_count += 1
            coverage.add(loc[0])

    return cover_count, coverage

def initialize_timestep_group(pf, reachable, loc_assignments, location_directory, pid_weight_map, radius):
    
    if len(reachable) == 0:
        pf_count = 0
        coverage = set()
    else:
        pf_count, coverage = cover_pid_group(pf, loc_assignments, location_directory, pid_weight_map, radius)
    
    return (pf, pf_count, coverage)

def cover_fac_sequential_r_group(potential_facilities, location_directory, pid_assignment, pid_weight_map, reachable_dict, k_facs, radius):
    
    paths = []
    total_coverage = 0
    
    START_OFFSET = min(pid_assignment.keys())
    END_TIME = max(pid_assignment.keys())
    
    dp = {START_OFFSET-1: {pf: {"count": 0, "path": []} for pf in potential_facilities}}
    fac_client = {}

    for k in range(0, k_facs):
        
        for hr in range(START_OFFSET, END_TIME+1):
            
            dp[hr] = {}
            loc_assignments = pid_assignment[hr]
            
            if k==0:
                results = [initialize_timestep_group(pf, [r for r in reachable_dict[pf] if r in location_directory.keys()], loc_assignments, location_directory, pid_weight_map, radius) for pf in potential_facilities]
                
                fac_client[hr] = {fac: {"count": count, "coverage": coverage} for fac, count, coverage in results}

            chosen_max = [fill_timestep(pf, [r for r in reachable_dict[pf] if r in location_directory.keys()], loc_assignments, dp[hr-1]) for pf in potential_facilities]
            
            for pf, max_count, max_path in chosen_max:
                dp[hr][pf] = {"count":max_count + fac_client[hr][pf]["count"],
                                 "path": max_path + [pf]}

        best_count, best_path = max([(dp[hr][fac]["count"], dp[hr][fac]["path"]) for fac in potential_facilities])
        paths.append(best_path)

        for i, fac in enumerate(best_path):
            
            h = i + START_OFFSET
            
            covered_clients = fac_client[h][fac]["coverage"]
            
            for pf in potential_facilities:
                fac_client[h][pf]["coverage"] = fac_client[h][pf]["coverage"].difference(covered_clients)
                fac_client[h][pf]["count"] = len(fac_client[h][pf]["coverage"])
            
            fac_client[h][fac]["coverage"] = set()
            fac_client[h][fac]["count"] = 0
        
        total_coverage += best_count
    
    return total_coverage, paths

def cover_fac_kary_group(potential_facilities, location_directory, pid_assignment, pid_weight_map, m, k_facs, lower_limit = 0.5, upper_limit = 10, n_jobs = 40):

    l = lower_limit
    r = upper_limit
    
    epsilon = 2e-3
    alpha = 1
    
    radius_paths = []
    opt_radius = r
    total_clients = sum(len(val) for val in pid_assignment.values())
    print("TOTAL CLIENTS: ", total_clients)
    
    reachable_dict = calculate_reachable(potential_facilities, location_directory, m)
    print("---CALCULATED REACHABLE FACILITIES---")
    
    while r-l > epsilon:
        
        partition = [l + i*(r - l)/(n_jobs+1) for i in range(n_jobs)]
        
        result = Parallel(n_jobs=n_jobs)(delayed(cover_fac_sequential_r_group)(potential_facilities, location_directory, pid_assignment, pid_weight_map, reachable_dict, k_facs, mid) for mid in partition)
        
        coverage = [result[i][0] for i in range(n_jobs)]
        
        if max(coverage) < total_clients:
            l = l + n_jobs*(r - l)/(n_jobs+1)
        else:
            for i in range(n_jobs):
                if coverage[i] >= total_clients:
                    radius_paths = result[i][1]
                    opt_radius = l + i*(r - l)/(n_jobs+1)
                    l = max(l, opt_radius - (r-l)/(n_jobs+1))
                    r = opt_radius
                    print(opt_radius)
                    break
    
    return opt_radius, radius_paths

def calculate_assignment_group(facilities, location_directory, pid_assignment, pid_weight_map):
    max_dist = 0
    
    group_assignment_costs = {}
    for val in pid_weight_map.values():
        group_assignment_costs[val] = {}
    
    START = min(pid_assignment.keys())
    END = max(pid_assignment.keys())
    
    for hr in range(START, END):
        for pid, client in pid_assignment[hr]:
            
            group = pid_weight_map[pid]
            
            coord_client = (location_directory[client]['latitude'], location_directory[client]['longitude'])
            client_cost = min([calculate_dist(coord_client, (location_directory[fac[hr-START]]['latitude'], location_directory[fac[hr-START]]['longitude'])) for fac in facilities])
            
            group_assignment_costs[group][pid] = client_cost
            
            if client_cost > max_dist:
                max_dist = client_cost
            
    return max_dist, group_assignment_costs

cover_runs = [run for run in os.listdir(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover") if run.endswith('.json') and "albe" in run]
random_seeds = set()

for file in cover_runs:
    with open(PROJECT_ROOT/"Experiments"/"output"/"budget_sensitivity"/"cover"/file) as f:
        seed = json.load(f)[0]
        random_seeds.add(seed)

for seed in random_seeds:
    
    county_name = "albemarle"
        
    potential_facilities, location_directory, pid_assignment = get_data(county_name, day, start, end, random_state = seed)
    print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))
    
    pid_weight_map = generate_weight_mapping(county_name, {(0, 10275):1, (10275, 41775):2, (41775, 1e10):3})
    
    print(seed)
    
    radius, facs = cover_fac_kary_group(potential_facilities, location_directory, pid_assignment, pid_weight_map, m, k_facs)
    obj, group_obj = calculate_assignment_group(facs, location_directory, pid_assignment, pid_weight_map)
    print("Weighted: ", obj)
    
    radius, facs = cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs)
    obj, unweighted_group_obj = calculate_assignment_group(facs, location_directory, pid_assignment, pid_weight_map)
    print("Unweighted: ", obj)
    
    with open(PROJECT_ROOT/"Experiments"/"output"/"weighted_inequality"/f"albe_weighted_cover_{seed}.json", "w") as f:
        json.dump({"seed": seed, "weighted": group_obj, "unweighted": unweighted_group_obj}, f)
