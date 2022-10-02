import geopy.distance
from itertools import chain, combinations
import tqdm
from joblib import Parallel, delayed

def calculate_dist(coord1, coord2):
    return geopy.distance.great_circle(coord1, coord2).km

def calculate_assignment_cost(facilities, location_directory, pid_assignment):
    max_dist = 0
    
    START = min(pid_assignment.keys())
    END = max(pid_assignment.keys())
    
    for hr in range(START, END):
        for pid, client in pid_assignment[hr]:
            coord_client = (location_directory[client]['latitude'], location_directory[client]['longitude'])
            client_cost = min([calculate_dist(coord_client, (location_directory[fac[hr-START]]['latitude'], location_directory[fac[hr-START]]['longitude'])) for fac in facilities])
            if client_cost > max_dist:
                max_dist = client_cost
    return max_dist