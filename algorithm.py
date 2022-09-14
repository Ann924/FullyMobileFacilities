import data
import tqdm
from joblib import Parallel, delayed
import geopy.distance

def calculate_dist(coord1, coord2):
    return geopy.distance.great_circle(coord1, coord2).km

def cover_pid_2(fac_loc, loc_assignment, location_directory, radius):
    
    cover_count = 0
    fac_coord = (location_directory[fac_loc]["latitude"], location_directory[fac_loc]["longitude"])
    
    coverage = set()
    
    for ind in range(len(loc_assignment)):
        loc = loc_assignment[ind]
        lid = loc[1]
        
        lid_coord = (location_directory[lid]["latitude"], location_directory[lid]["longitude"])
        if calculate_dist(fac_coord, lid_coord) <= radius:
            cover_count += 1
            coverage.add(loc[0])

    return cover_count, coverage

def initialize_timestep_2(pf, reachable, loc_assignments, location_directory, radius):
    
    if len(reachable) == 0:
        pf_count = 0
        coverage = set()
    else:
        pf_count, coverage = cover_pid_2(pf, loc_assignments, location_directory, radius)
    
    return (pf, pf_count, coverage)

def fill_timestep_2(pf, reachable, loc_assignments, dp):
    
    if len(reachable) == 0:
        max_count = 0
        max_path = set()
    else:
        max_count, max_path = max((dp[adj]["count"], dp[adj]["path"]) for adj in reachable)
    
    return (pf, max_count, max_path)

def cover_fac_parallel_r_2(potential_facilities, location_directory, pid_assignment, reachable_dict, k_facs, radius):
    
    paths = []
    total_coverage = 0
    
    START_OFFSET = min(pid_assignment.keys())
    END_TIME = max(pid_assignment.keys())
    
    dp = {START_OFFSET-1: {pf: {"count": 0, "path": []} for pf in potential_facilities}}
    fac_client = {}

    for k in tqdm.tqdm(range(0, k_facs)):
        
        for hr in range(START_OFFSET, END_TIME+1):
            
            dp[hr] = {}
            loc_assignments = pid_assignment[hr]
            
            if k==0:
                results = Parallel(n_jobs=5)(delayed(initialize_timestep_2)
                                            (pf, [r for r in reachable_dict[pf] if r in location_directory.keys()],
                                             loc_assignments, location_directory, radius)
                                                 for pf in potential_facilities)
                fac_client[hr] = {fac: {"count": count, "coverage": coverage} for fac, count, coverage in results}

            chosen_max = Parallel(n_jobs=5)(delayed(fill_timestep_2)(pf, [r for r in reachable_dict[pf] if r in location_directory.keys()], loc_assignments, dp[hr-1])
                                            for pf in potential_facilities)
            
            for pf, max_count, max_path in chosen_max:
                dp[hr][pf] = {"count":max_count + fac_client[hr][pf]["count"],
                                 "path": max_path + [pf]}

        best_count, best_path = max([(dp[hr][fac]["count"],
                                      dp[hr][fac]["path"]) for fac in potential_facilities])
        paths.append(best_path)

        for i, fac in enumerate(best_path):
            
            h = i + START_OFFSET
            
            covered_clients = fac_client[h][fac]["coverage"]
            #print(fac, len(covered_clients))
            
            for pf in potential_facilities:
                fac_client[h][pf]["coverage"] = fac_client[h][pf]["coverage"].difference(covered_clients)
                fac_client[h][pf]["count"] = len(fac_client[h][pf]["coverage"])
            
            fac_client[h][fac]["coverage"] = set()
            fac_client[h][fac]["count"] = 0
        
        print(best_count, best_path)
        
        total_coverage += best_count
    
    return total_coverage, paths

def reachable_helper(potential_facilities, location_directory, lid_1, m):
    reachable = []
    
    for lid_2 in potential_facilities:
        
        coord_1 = (location_directory[lid_1]["latitude"], location_directory[lid_1]["longitude"])
        coord_2 = (location_directory[lid_2]["latitude"], location_directory[lid_2]["longitude"])

        dist = calculate_dist(coord_1, coord_2)
        if dist<=m:
            reachable.append(lid_2)
    
    return (lid_1, reachable)

def calculate_reachable(potential_facilities, location_directory, m):
    
    
    results = Parallel(n_jobs=5)(delayed(reachable_helper)
                            (potential_facilities, location_directory, pf, m)
                                 for pf in potential_facilities)
    
    reachable_dict = {pf: reachable for pf,reachable in results}
    return reachable_dict
        

def cover_fac_binary_2(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.5, upper_limit=3):

    l = lower_limit
    r = upper_limit
    
    epsilon = 1e-3
    alpha = 1
    
    radius_paths = []
    opt_radius = r
    total_clients = sum(len(val) for val in pid_assignment.values())
    print("TOTAL CLIENTS: ", total_clients)
    
    reachable_dict = calculate_reachable(potential_facilities, location_directory, m)
    print("---CALCULATED REACHABLE FACILITIES---")
    
    while r-l > 1e-3:
        
        mid = (l+r)/2
        
        client_coverage, paths = cover_fac_parallel_r_2(potential_facilities, location_directory,
                                                        pid_assignment, reachable_dict, k_facs, mid)
        print(mid, client_coverage)
        
        if client_coverage < total_clients:
            l = mid
        else:
            radius_paths = paths
            opt_radius = mid
            r = mid
    
    return opt_radius, radius_paths