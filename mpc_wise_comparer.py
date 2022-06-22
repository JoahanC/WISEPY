from visualizer_functions import visualizer
from mpc_wise_functions import *
import sys

new_epochs = comparer(sys.argv[1], sys.argv[2], False)
visualizer(sys.argv[1], sys.argv[2], sys.argv[3])

print("GATOR source id query string:")
sid_string = "source_id in ("
counter = 0
for epoch in new_epochs:
    sid_string += "'" + new_epochs[epoch][0] + "', "
print(sid_string[:-2] + ')')