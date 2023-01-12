import math
import numpy as np
from ortools.linear_solver import pywraplp
from ortools.linear_solver.pywraplp import Constraint, Variable, Objective
from FullyMobile.utils import calculate_dist
import time

class FlowLP():
    def __init__(self, potential_facilities, location_directory, pid_assignment, m, k_facs, radius, solver_id: str = "GLOP"):

        self.solver = pywraplp.Solver.CreateSolver(solver_id)
        
        self.potential_facilities = potential_facilities
        self.location_directory = location_directory
        self.pid_assignment = pid_assignment
        self.k = k_facs
        self.m = m
        self.radius = radius
        self.max_timestep = max(pid_assignment.keys())

        self.client_balls = self.calculate_client_balls()

        # Partial evaluation storage
        self.partials = {}

        start_var = time.time()
        self.init_variables()
        end_var = time.time()
        
        print(end_var-start_var)
        
        start_var = time.time()
        self.init_constraints()
        end_var = time.time()
        
        print(end_var-start_var)

        
    def calculate_client_balls(self, n_jobs = 40):
        
        client_balls = {}
        
        for time in self.pid_assignment.keys():
            client_balls[time] = {}
            for pid, pid_loc in self.pid_assignment[time]:
                client_balls[time][pid_loc] = []
                
                for pf in self.potential_facilities:
                    coord_1 = (self.location_directory[pf]["latitude"], self.location_directory[pf]["longitude"])
                    coord_2 = (self.location_directory[pid_loc]["latitude"], self.location_directory[pid_loc]["longitude"])
                    
                    dist = calculate_dist(coord_1, coord_2)
                    if dist<=self.radius:
                        client_balls[time][pid_loc].append(pf)
        
        return client_balls
                        
    def init_variables(self):
        """Declare variables as needed"""
        
        # indicator variable for whether facility is in solution
        self.Y = {}
        
        # indicator variable for whether client is serviced by a facility
        self.X = {}
        
        # indicator variable for movement between timesteps
        self.Z = {}
        
        for time in self.pid_assignment.keys():
            self.Y[time] = {}
            self.X[time] = {}
            self.Z[time] = {}
            
            for fac in self.potential_facilities:
                self.Y[time][fac] = self.solver.NumVar(0, 1, f"y_{time}{fac}")
                self.X[time][fac] = {}
                self.Z[time][fac] = {}
                
                for pid, pid_loc in self.pid_assignment[time]:
                    if pid_loc not in self.X[time][fac].keys():
                        self.X[time][fac][pid_loc] = self.solver.NumVar(0, 1, f"x_{time}{fac}{pid_loc}")
                
                if time!=self.max_timestep:
                    for next_fac in self.potential_facilities:
                        self.Z[time][fac][next_fac] = self.solver.NumVar(0, 1, f"z_{time}{fac}{next_fac}")            

    def init_constraints(self):
        """Initializes the constraints according to the relaxed LP formulation of MinExposed"""

        # budget constraint for each timestep
        self.budget_constraints = {}
        
        # constraint for facility and client assignment
        self.client_constraints = {}
        
        # constraint for movement
        self.leave_constraints = {}
        self.enter_constraints = {}

        for time in self.Y.keys():
            self.budget_constraints[time] = self.solver.Constraint(0, self.k, f"budget_{time}")
            for fac, y_var in self.Y[time].items():
                self.budget_constraints[time].SetCoefficient(y_var, 1)
                
                # constraint for facility assignment
                for pid_loc, x_var in self.X[time][fac].items():
                    self.solver.Add(x_var <= y_var)
                
            self.client_constraints[time] = {}
            for pid, pid_loc in self.pid_assignment[time]:
                if pid_loc not in self.client_constraints[time].keys():
                    self.client_constraints[time][pid_loc] = self.solver.Constraint(1, 1, f"client_assignment_{time}_{pid_loc}")
                
                    for pf in self.client_balls[time][pid_loc]:
                        self.client_constraints[time][pid_loc].SetCoefficient(self.X[time][pf][pid_loc], 1)

            if time != self.max_timestep:
                
                self.leave_constraints[time] = {}
                self.enter_constraints[time] = {}
                
                for fac, y_var in self.Y[time].items():
                    #self.leave_constraints[fac] = self.solver.Constraint(y_var, y_var, f"leave_constraint_{time}_{fac}")
                    self.leave_constraints[fac] = self.solver.Add(self.solver.Sum([z_var for z_var in self.Z[time][fac].values()]) == y_var)
                    for next_fac, z_var in self.Z[time][fac].items():
                        #self.leave_constraints[fac].SetCoefficient(z_var, 1)
                        
                        coord_1 = (self.location_directory[fac]["latitude"], self.location_directory[fac]["longitude"])
                        coord_2 = (self.location_directory[next_fac]["latitude"], self.location_directory[next_fac]["longitude"])
    
                        dist = calculate_dist(coord_1, coord_2)
                        if dist > self.m:
                            self.solver.Add(self.Z[time][fac][next_fac] == 0)
                
                for next_fac, y_var in self.Y[time+1].items():
                    self.enter_constraints[next_fac] = self.solver.Add(self.solver.Sum([self.Z[time][fac][next_fac] for fac in self.Z[time].keys()]) == y_var)
                    #self.solver.Constraint(y_var, y_var, f"enter_constraint_{time}_{next_fac}")
                    #for fac in self.Z[time].keys():
                    #    self.enter_constraints[next_fac].SetCoefficient(self.Z[time][fac][next_fac], 1)

    def solve_lp(self):
        
        self.status = self.solver.Solve()
        
        if self.status == pywraplp.Solver.OPTIMAL:
            print('Problem solved in %f milliseconds' % self.solver.wall_time())
            print('Problem solved in %d iterations' % self.solver.iterations())
        else:
            print('The problem does not have an optimal solution.')

    def get_variable_solution(self):
        if self.status == pywraplp.Solver.OPTIMAL:
            y = {time: {fac: self.Y[time][fac].solution_value() for fac in self.Y[time].keys()} for time in self.Y.keys()}
            x = {time: {fac: {pid: self.X[time][fac][pid].solution_value() for pid in self.X[time][fac].keys()} for fac in self.X[time].keys()} for time in self.X.keys()}
            z = {time: {fac: {next_fac: self.Z[time][fac][next_fac].solution_value() for next_fac in self.Z[time][fac].keys()} for fac in self.Z[time].keys()} for time in self.Z.keys()}
            return y, x, z
        else:
            print("Not optimal")
            return {}, {}, {}
    
    