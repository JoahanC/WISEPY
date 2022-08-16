from asteroid import Asteroid
from ds9_interface_functions import generate_script
import os

files = os.listdir("./database_files/wise")
query = {}
for file in files:
    packed_name = file[:5]
    band = file[6]
    if packed_name not in query:
        query[packed_name] = [int(band)]
    else:
        query[packed_name].append(int(band))

temp = query["G1989"]
del query["G1989"]
query[161989] = temp

formatted_query = {}
for name in query:
    formatted_query[int(name)] = query[name]

for name in formatted_query:
    print(f"Running WISEPY for {name}.")
    asteroid = Asteroid(name, formatted_query[name])
    asteroid.run_comparison()

"""cacus = Asteroid(161989, [2, 4])
#print(cacus.observations)
#print(cacus.ra_dec_values)
#print(cacus.band_images)
#print(cacus.unique_instances)
#cacus.display_unique_observations()
#cacus.generate_new_table(2)
cacus.generate_new_table(4)
cacus.generate_full_table(2)
#cacus.generate_script(2, 1, 10)
cacus.generate_flux_snr_plots()
#generate_script("161989", 2, 1, 10)
#print(cacus.split_wise_observations[2])
#cacus.load_files("G1989_2band.txt", 2)"""