from visualizer_functions import visualizer
from mpc_wise_functions import *
import sys

new_epochs = comparer(sys.argv[1], sys.argv[2], False)
visualizer(sys.argv[1], sys.argv[2], sys.argv[3])

print("All new Source Ids:")
sid_string = ""
counter = 0
for epoch in new_epochs:
    print(new_epochs[epoch])
    counter += 1
    if counter <= 50:
        sid_string += "'" + new_epochs[epoch][0] + "', "
print(sid_string[:-2])
