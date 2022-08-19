"""
This file runs WISEPY on all valid inputted asteroids.
"""
import os
from asteroid import Asteroid

# Runs wisepy on all current objects.
files = os.listdir("./database_files/wise")
query = {}
for file in files:
    packed_name = file[:5]
    band = file[6]
    if packed_name not in query:
        query[packed_name] = [int(band)]
    else:
        query[packed_name].append(int(band))

formatted_query = {}
for name in query:
    formatted_query[int(name)] = query[name]

for name in formatted_query:
    print(f"Running WISEPY for {name}.")
    asteroid = Asteroid(name, formatted_query[name])
    asteroid.run_comparison()