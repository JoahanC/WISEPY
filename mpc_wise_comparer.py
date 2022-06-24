from visualizer_functions import visualizer
from mpc_wise_functions import *
import sys


mpc_file = "input_data/" + sys.argv[1]
wise_file = "input_data/" + sys.argv[2]

new_epochs = comparer(mpc_file, wise_file, False)
visualizer(mpc_file, wise_file, sys.argv[3])

print("GATOR source id query string:")
sid_string = "source_id in ("
counter = 0
for epoch in new_epochs:
    sid_string += "'" + new_epochs[epoch][0] + "', "
print(sid_string[:-2] + ')')