import tqdm
from joblib import Parallel, delayed
from FullyMobile.utils import calculate_dist
import geopy.distance


"""
K-Supplier Baseline
"""
def k_supplier(potential_facilities, location_directory, pid_assignment, k: int, lower_limit = 0.1, upper_limit = 20):
    """
    Solves k-supplier (where client locations and facility locations may not overlap) with Hochbaum-Shmoys
    3-approximation algorithm
    """
    
    clients = set()
    for hr, vals in pid_assignment.items():
        for val in vals:
            clients.add(val[1])
    
    clients = list(clients)
    
    l = lower_limit
    r = upper_limit

    to_ret = -1
    EPSILON = 2e-3
    
    while r-l > EPSILON:
    
        mid = l + (r - l) / 2

        if len(_check_radius(location_directory, mid, clients)) <= k:
            facilities: List[int] = _locate_facilities(location_directory, mid,
                                    _check_radius(location_directory, mid, clients), potential_facilities, k)
            if facilities:
                to_ret = mid
                r = mid
            else:
                l = mid
        else:
            l = mid
    
    return _locate_facilities(location_directory, to_ret, _check_radius(location_directory, to_ret, clients), potential_facilities, k)

def _check_radius(location_directory, radius: int, clients):
    """Determine the maximal independent set of pairiwse independent client balls with given radius
    
    RETURNS
    ----------
    pairwise_disjoint
        maximal independent pairwise disjoint set of clients, where disjoint is defined as greater than a distance
        of 2*radius apart
    """
    
    pairwise_disjoint = set()

    V = set(clients)
    while len(V)!=0:
        v = V.pop()
        pairwise_disjoint.add(v)
        
        remove = set()
        for i in V:
            coord1 = (location_directory[v]["latitude"], location_directory[v]["longitude"])
            coord2 = (location_directory[i]["latitude"], location_directory[i]["longitude"])
            if calculate_dist(coord1, coord2) <= 2*radius:
                remove.add(i)
        V-=remove
    
    return pairwise_disjoint

def _locate_facilities(location_directory, radius: int, pairwise_disjoint, locations, k: int):
    """Select a facility to open within the given radius for each pairwise_disjoint client
    """
    
    facilities = set()
    for c in pairwise_disjoint:
        for l in locations:
            coord1 = (location_directory[c]["latitude"], location_directory[c]["longitude"])
            coord2 = (location_directory[l]["latitude"], location_directory[l]["longitude"])
            if calculate_dist(coord1, coord2) <= radius:
                facilities.add(l)
                break
    
    if len(facilities) < len(pairwise_disjoint):
        return None
    
    #Check if k larger than the number of possible facility locations
    k = min(k, len(locations))
    
    #Randomly add facilities for leftover budget
    if k>len(facilities):
        unopened_facilities = set(locations)-facilities
        for i in range(k-len(facilities)):
            facilities.add(unopened_facilities.pop())
    
    return list(facilities)

def stationary_supplier(potential_facilities, location_directory, pid_assignment, k_facs):
    clients = set()
    for key, vals in pid_assignment.items():
        clients = clients.union(vals)
    
    facs = k_supplier(potential_facilities, location_directory, pid_assignment, k_facs)
    return [[f]*len(pid_assignment.keys()) for f in facs]

"""
Cover Algorithm
"""

def cover_pid(fac_loc, loc_assignment, location_directory, radius):
    
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

def initialize_timestep(pf, reachable, loc_assignments, location_directory, radius):
    
    if len(reachable) == 0:
        pf_count = 0
        coverage = set()
    else:
        pf_count, coverage = cover_pid(pf, loc_assignments, location_directory, radius)
    
    return (pf, pf_count, coverage)

def fill_timestep(pf, reachable, loc_assignments, dp):
    
    if len(reachable) == 0:
        max_count = 0
        max_path = set()
    else:
        max_count, max_path = max((dp[adj]["count"], dp[adj]["path"]) for adj in reachable)
    
    return (pf, max_count, max_path)

def cover_fac_sequential_r(potential_facilities, location_directory, pid_assignment, reachable_dict, k_facs, radius):
    
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
            
            # TODO: simplify and refactor this process #
            if k==0:
                results = [initialize_timestep(pf, [r for r in reachable_dict[pf] if r in location_directory.keys()], loc_assignments, location_directory, radius) for pf in potential_facilities]
                
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

def reachable_helper(potential_facilities, location_directory, lid_1, m):
    reachable = []
    
    for lid_2 in potential_facilities:
        
        coord_1 = (location_directory[lid_1]["latitude"], location_directory[lid_1]["longitude"])
        coord_2 = (location_directory[lid_2]["latitude"], location_directory[lid_2]["longitude"])

        dist = calculate_dist(coord_1, coord_2)
        if dist<=m:
            reachable.append(lid_2)
    
    return (lid_1, reachable)

def calculate_reachable(potential_facilities, location_directory, m, n_jobs = 40):
    results = Parallel(n_jobs=n_jobs)(delayed(reachable_helper)
                            (potential_facilities, location_directory, pf, m)
                                 for pf in potential_facilities)
    
    reachable_dict = {pf: reachable for pf,reachable in results}
    return reachable_dict

def cover_fac_kary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.5, upper_limit = 3, n_jobs = 40):

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
        
        result = Parallel(n_jobs=n_jobs)(delayed(cover_fac_sequential_r)(potential_facilities, location_directory, pid_assignment, reachable_dict, k_facs, mid) for mid in partition)
        
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
    
########################## UNUSED ALGORITHMS #################################

### Binary Search ###
"""
def cover_fac_parallel_r(potential_facilities, location_directory, pid_assignment, reachable_dict, k_facs, radius, n_jobs=40):
    
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
                results = Parallel(n_jobs=n_jobs)(delayed(initialize_timestep)
                                            (pf, [r for r in reachable_dict[pf] if r in location_directory.keys()],
                                             loc_assignments, location_directory, radius)
                                                 for pf in potential_facilities)
                fac_client[hr] = {fac: {"count": count, "coverage": coverage} for fac, count, coverage in results}

            chosen_max = Parallel(n_jobs=n_jobs)(delayed(fill_timestep)(pf, [r for r in reachable_dict[pf] if r in location_directory.keys()], loc_assignments, dp[hr-1])
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
"""

"""

def cover_fac_binary(potential_facilities, location_directory, pid_assignment, m, k_facs, lower_limit = 0.5, upper_limit = 3, n_jobs=40):

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
        
        mid = (l+r)/2
        
        client_coverage, paths = cover_fac_parallel_r(potential_facilities, location_directory,
                                                        pid_assignment, reachable_dict, k_facs, mid, n_jobs=n_jobs)
        print(mid, client_coverage)
        
        if client_coverage < total_clients:
            l = mid
        else:
            radius_paths = paths
            opt_radius = mid
            r = mid
    
    return opt_radius, radius_paths
"""