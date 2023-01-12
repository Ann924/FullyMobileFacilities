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
old_assign = pid_assignment
print(len(potential_facilities), len(location_directory), sum(len(val) for val in pid_assignment.values()))

from FullyMobile.lp import *
lp = FlowLP(potential_facilities, location_directory, pid_assignment, 5, 5, 5)

start = time.time()
lp.solve_lp()
end = time.time()

with open("lp_run.json", 'w') as f:
    json.dump(end-start, f)