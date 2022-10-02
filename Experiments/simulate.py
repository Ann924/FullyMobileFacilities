import data
import algorithm
import time

def run():
    
    potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", 5, 6, 20)

    k_count = 5
    m = 6
    opt_radius, paths = algorithm.cover_fac_binary_2(potential_facilities, location_directory, pid_assignment, m, k_count, lower_limit = 0.5, upper_limit=3)
    
    for k in range(k_count):
        print(f"Facility #{k}: ", str(paths[k]))
    
    print(f"Final Radius: {opt_radius}")

start = time.time()
run()
end = time.time()
print("Total Time: ", end-start)