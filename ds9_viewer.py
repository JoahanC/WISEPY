from mpc_wise_functions import *
from ds9_interface_functions import *
import os
import re
import sys


mpc_file = "input_data/" + sys.argv[1] + ".txt"
wise_file = "input_data/" + sys.argv[1] + ".tbl"

source_ids = generate_source_ids_list(mpc_file, wise_file)
view_all = input(f"There are {len(source_ids) + 1} unique epochs detected.\n"
                 + "Would you like to look at all of them? [Y]/[N]\n")
if view_all.capitalize() == 'Y':
    script = generate_script(source_ids, 0, len(source_ids) + 1, sys.argv[1])
    os.popen(script)
else:
    lower_bound = input(f"Which image out of 1-{len(source_ids) + 1} would "
                        + "you want to start at\n")
    upper_bound = input(f"Which image out of 1-{len(source_ids) + 1} would "
                        + "you want to end at\n")
    script = generate_script(source_ids, int(lower_bound) - 1, 
                            int(upper_bound), sys.argv[1])
    os.popen(script)