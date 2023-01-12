import tqdm
from joblib import Parallel, delayed
from FullyMobile.utils import calculate_dist
import geopy.distance
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.sparse.csgraph import maximum_flow
import numpy as np

"""
K-Supplier Iterative Baseline
"""

def check_matching(potential_facilities, location_directory, previous_facilities, timestep_balls, k, m, radius):
    if len(previous_facilities) == 0:
        facilities: List[int] = _locate_facilities(location_directory, radius, timestep_balls, potential_facilities, k)
        return facilities
    else:
        distance_dict = {}
        
        data = []
        indices = []
        indptr = [0]
        
        for prev_iter, prev in enumerate(previous_facilities):
            coord_prev = (location_directory[prev]["latitude"], location_directory[prev]["longitude"])
            
            edge_count = 0
            
            for ball_iter, next_ball in enumerate(timestep_balls):
                
                coord_client = (location_directory[next_ball]["latitude"], location_directory[next_ball]["longitude"])
                
                for fac in potential_facilities:
                    coord_fac = (location_directory[fac]["latitude"], location_directory[fac]["longitude"])
                    if calculate_dist(coord_client, coord_fac) <= radius and calculate_dist(coord_prev, coord_fac) <= m:
                        distance_dict[(prev_iter, ball_iter)] = fac
                        
                        edge_count += 1
                        
                        data.append(1)
                        indices.append(ball_iter)
                        
                        break
            
            indptr.append(indptr[-1] + edge_count)
        
        #print(data, indices, indptr)
        
        graph = csr_matrix((data, indices, indptr), shape=(len(previous_facilities), len(timestep_balls)))
        match = maximum_bipartite_matching(graph, perm_type = 'column')
        
        if len(timestep_balls) >= len(previous_facilities) and (match == -1).sum() > 0:
            return {}
        else:
            fac = []
            for row, col in enumerate(match):
                if col!=-1:
                    fac.append(distance_dict[(row, col)])
                else:
                    fac.append(previous_facilities[row])
            return fac
        
def k_supplier_step(potential_facilities, location_directory, k, m, previous_facilities, clients, lower_limit = 0.1, upper_limit = 20):
    
    l = lower_limit
    r = upper_limit

    to_ret = -1
    ret_facilities = []
    EPSILON = 2e-3
    
    while r-l > EPSILON:
    
        mid = l + (r - l) / 2
        
        pairwise_disjoint = _check_radius(location_directory, mid, clients)

        if len(pairwise_disjoint) <= k:
            facilities = check_matching(potential_facilities, location_directory, previous_facilities, pairwise_disjoint, k, m, mid)
            
            if facilities:
                to_ret = mid
                r = mid
                ret_facilities = facilities
            else:
                l = mid
        
        else:
            l = mid
    
    return to_ret, ret_facilities

def iterative_supplier(potential_facilities, location_directory, pid_assignment, k: int, m: float, lower_limit = 0.1, upper_limit = 20):
    
    previous_facilities = []
    facility_paths = [[] for _ in range(k)]
    
    for hr, vals in pid_assignment.items():
        
        clients = set()
        for val in vals:
            clients.add(val[1])
        
        radius, timestep_facilities = k_supplier_step(potential_facilities, location_directory, k, m, previous_facilities, clients, lower_limit = lower_limit, upper_limit = upper_limit)
        print(radius, timestep_facilities)
        for i in range(len(facility_paths)):
            facility_paths[i].append(timestep_facilities[i])
        previous_facilities = timestep_facilities
    
    return facility_paths

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

def stationary_supplier(potential_facilities, location_directory, pid_assignment, k_facs, lower_limit = 0.1, upper_limit = 20):

    facs = k_supplier(potential_facilities, location_directory, pid_assignment, k_facs, lower_limit = lower_limit, upper_limit = upper_limit)
    
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


"""
Flow Algorithm
"""

def flow_search(potential_facilities, location_directory, pid_assignment, k, m, lower_limit = 0.1, upper_limit = 20):
        
    l = lower_limit
    r = upper_limit

    to_ret = -1
    ret_facilities = []
    EPSILON = 2e-3
    
    reachable = calculate_reachable(potential_facilities, location_directory, m)
    
    while r-l > EPSILON:
    
        mid = l + (r - l) / 2
        
        client_balls = {}
        good_radius = True
        
        for t in pid_assignment.keys():
            clients = {c for pid,c in pid_assignment[t]}
            pairwise_disjoint = _check_radius(location_directory, mid, clients)
            
            if len(pairwise_disjoint)>k:
                good_radius = False
                break
            
            client_balls[t] = pairwise_disjoint

        if good_radius:
            facilities = build_flow(potential_facilities, location_directory, client_balls, reachable, k, m, mid)

            if facilities:
                to_ret = mid
                r = mid
                ret_facilities = facilities
                #print(mid, len(facilities))
            else:
                #print(mid, "INFEASIBLE")
                l = mid
        
        else:
            #print(mid, "IND SETS")
            l = mid
    
    return ret_facilities

def build_flow(potential_facilities, location_directory, client_balls, reachable, k, m, radius):
    
    fac_to_clients = {}
    
    for pf in potential_facilities:
        fac_to_clients[pf] = {}
        coord1 = (location_directory[pf]["latitude"], location_directory[pf]["longitude"])
        
        for t, clients in client_balls.items():
            fac_to_clients[pf][t] = set()
            
            for c in clients:
                coord2 = (location_directory[c]["latitude"], location_directory[c]["longitude"])
                dist = calculate_dist(coord1, coord2)
                if dist <= radius:
                    fac_to_clients[pf][t].add(c)
    
    #print('calculated reachable facilities')

    data = []
    indices = []
    indptr = [0]
    
    node_to_ind = {}
    count = 0
    for t in client_balls.keys():
        node_to_ind[t] = {"fac_in": {}, "client_in": {}, "client_out": {},"fac_out":{}}
        for pf in potential_facilities:
            node_to_ind[t]["fac_in"][pf] = count
            count += 1
        for c in client_balls[t]:
            node_to_ind[t]["client_in"][c] = count
            count += 1
        for c in client_balls[t]:
            node_to_ind[t]["client_out"][c] = count
            count += 1
        for pf in potential_facilities:
            node_to_ind[t]["fac_out"][pf] = count
            count += 1
    
    #s, t, s', t' have not been accounted for yet...
    S_STAR = count
    S_NODE = S_STAR+1
    T_NODE = S_NODE+1
    S_PRIME = T_NODE+1
    T_PRIME = S_PRIME+1
    
    row_count = 0
    for t in client_balls.keys():
        client_to_fac = {}
        for pf in potential_facilities:
            row_count += 1
            connection_count = 0
            
            if len(fac_to_clients[pf][t])== 0:
                
                connection_count += 1
                data.append(k)
                indices.append(node_to_ind[t]["fac_out"][pf])
                
            for c in fac_to_clients[pf][t]:
                
                if c not in client_to_fac.keys():
                    client_to_fac[c] = [pf]
                else:
                    client_to_fac[c].append(pf)
                
                connection_count += 1
                indices.append(node_to_ind[t]["client_in"][c])
                data.append(k)
            
            indptr.append(indptr[-1] + connection_count)
        
        for c in client_balls[t]:
            row_count += 1
            # may be another larger capacity...
            data.append(k)
            indices.append(node_to_ind[t]["client_out"][c])
            
            # add connection to t'
            data.append(1)
            indices.append(T_PRIME)
            indptr.append(indptr[-1]+2)
        
        for c in client_balls[t]:
            row_count += 1
            connection_count = 0
            for fac in client_to_fac[c]:
                indices.append(node_to_ind[t]["fac_out"][fac])
                data.append(k)
                connection_count += 1
            indptr.append(indptr[-1]+connection_count)
        
        for pf in potential_facilities:
            row_count += 1
            if t == max(client_balls.keys()):
                #add connections to t
                data.append(k)
                indices.append(T_NODE)
                indptr.append(indptr[-1]+1)
            else:
                connection_count = 0
                for fac in reachable[pf]:
                    indices.append(node_to_ind[t+1]["fac_in"][fac])
                    data.append(k)
                    connection_count += 1
                indptr.append(indptr[-1]+connection_count)
    
    saturation = 0
    
    # add constraint on total flow
    data.append(k)
    indices.append(T_PRIME)
    indptr.append(indptr[-1]+1)
    
    # add s edges
    row_count += 1
    data += [k for _ in range(len(potential_facilities))]
    indices += [node_to_ind[min(client_balls.keys())]["fac_in"][fac] for fac in potential_facilities]
    indptr.append(indptr[-1]+len(potential_facilities))
    
    # add t edge (going back into s*)
    row_count += 1
    data.append(k)
    indices.append(S_STAR)
    indptr.append(indptr[-1]+1)
    
    # add s' edges
    row_count += 1
    connection_count = 0
    for t in client_balls.keys():
        for client, indx in node_to_ind[t]["client_out"].items():
            connection_count += 1
            indices.append(indx)
            data.append(1)
    
    data.append(k)
    indices.append(S_NODE)
    
    indptr.append(indptr[-1]+connection_count + 1)
    saturation = connection_count + k
    
    # add t' edges
    row_count += 1
    indptr.append(indptr[-1])
    
    graph = csr_matrix((data, indices, indptr), shape=(T_PRIME+1, T_PRIME+1))
    flow = maximum_flow(graph, S_PRIME, T_PRIME, method='edmonds_karp')
    if flow.flow_value == saturation:
        #print("FEASIBLE!")
        #fac = open_facilities(graph.toarray()[:,:S_STAR], flow.flow.toarray()[:,:S_STAR], node_to_ind, fac_to_clients, potential_facilities, client_balls, S_NODE)
        fac = open_facilities(graph[:,:S_STAR], flow.flow[:,:S_STAR], node_to_ind, fac_to_clients, potential_facilities, client_balls, S_NODE)
        return fac
    else:
        return


def open_facilities(graph, flow, node_to_ind, fac_to_clients, potential_facilities, client_balls, S_NODE):
    
    facility_paths = []
    
    while len((flow[S_NODE]>0).indices) > 0:
        
        path = []
        
        fac = (flow[S_NODE]>0).indices[0]
        flow[S_NODE, fac]-=1
        
        t = min(client_balls.keys())
        ind_to_node = {ind: node for node, ind in node_to_ind[t]["fac_in"].items()}
        
        path.append(ind_to_node[fac])
        
        while len(path)<len(client_balls.keys()):
            
            t += 1
            ind_to_node = {ind: node for node, ind in node_to_ind[t]["fac_in"].items()}
            
            previous_facility = path[-1]
            index = node_to_ind[t-1]["fac_in"][previous_facility]
            next_step = (flow[index]>0).indices[0]
            flow[index, next_step]-=1

            # does not flow through client ball
            if len(fac_to_clients[previous_facility][t-1]) == 0:

                next_fac = (flow[next_step]>0).indices[0]

                flow[next_step, next_fac]-=1

            # it flows through a client ball
            else:
                next_client_ball = (graph[next_step]>0).indices[0]
                exit_node = (flow[next_client_ball]>0).indices[0]
                next_fac = (flow[exit_node]>0).indices[0]

                flow[next_client_ball, exit_node]-=1
                flow[exit_node, next_fac]-=1
        
            path.append(ind_to_node[next_fac])
        
        facility_paths.append(path)
        #print(facility_paths)
        
    return facility_paths
    

'''def open_facilities(graph, flow, node_to_ind, fac_to_clients, potential_facilities, client_balls, S_NODE):
    
    facility_paths = []
    
    while len(np.where(flow[S_NODE]>0)[0]) > 0:
        
        path = []
        
        fac = np.where(flow[S_NODE]>0)[0][0]
        flow[S_NODE, fac]-=1
        
        t = min(client_balls.keys())
        ind_to_node = {ind: node for node, ind in node_to_ind[t]["fac_in"].items()}
        
        path.append(ind_to_node[fac])
        
        while len(path)<len(client_balls.keys()):
            
            t += 1
            ind_to_node = {ind: node for node, ind in node_to_ind[t]["fac_in"].items()}
            
            previous_facility = path[-1]
            index = node_to_ind[t-1]["fac_in"][previous_facility]
            next_step = np.where(flow[index]>0)[0][0]
            flow[index, next_step]-=1

            # does not flow through client ball
            if len(fac_to_clients[previous_facility][t-1]) == 0:

                next_fac = np.where(flow[next_step]>0)[0][0]

                flow[next_step, next_fac]-=1

            # it flows through a client ball
            else:
                next_client_ball = np.where(graph[next_step]>0)[0][0]
                exit_node = np.where(flow[next_client_ball]>0)[0][0]
                next_fac = np.where(flow[exit_node]>0)[0][0]

                flow[next_client_ball, exit_node]-=1
                flow[exit_node, next_fac]-=1
        
            path.append(ind_to_node[next_fac])
        
        facility_paths.append(path)
        
    return facility_paths'''

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