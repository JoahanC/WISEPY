from mpc_wise_functions import *
from ds9_interface_functions import *
import os
import sys

mpc_code = sys.argv[1]
bands = int(sys.argv[2])

source_ids = generate_source_ids_list(mpc_code, bands)

view_all = input(f"There are {len(source_ids)} unique epochs detected.\n"
                 + "Would you like to look at all of them? [Y]/[N]\n")

if view_all.capitalize() == 'Y':
    script = generate_script(source_ids, 0, len(source_ids) * int(sys.argv[2]) + 1, sys.argv[1], int(sys.argv[2]))
    os.popen(script)

else:
    lower_bound = input(f"Which image out of 1-{len(source_ids) * int(sys.argv[2])} would "
                        + "you want to start at\n")
    upper_bound = input(f"Which image out of 1-{len(source_ids) * int(sys.argv[2])} would "
                        + "you want to end at\n")
    script = generate_script(source_ids, int(lower_bound) - 1, 
                            int(upper_bound), sys.argv[1], int(sys.argv[2]))
    os.popen(script)