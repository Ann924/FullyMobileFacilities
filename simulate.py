import data
import algorithm

def run():
    potential_facilities, location_directory, pid_assignment = data.get_data("charlottesville_city", 5, 6, 20)

    k_count = 5
    m = 6
    opt_radius, paths = algorithm.cover_fac_binary_2(potential_facilities, location_directory, pid_assignment, m, k_count, lower_limit = 0.5, upper_limit=2)
    
    for k in range(k_count):
        print(f"Facility #{k}: ", str(paths[k]))
    
    print(f"Final Radius: {opt_radius}")

run()