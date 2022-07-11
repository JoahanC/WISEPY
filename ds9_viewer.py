from mpc_wise_functions import *
from ds9_interface_functions import *
import os
import re
import sys

if int(sys.argv[2]) == 2:
    mpc_file = "input_data/" + sys.argv[1] + ".txt"
    wise_file = "input_data/" + sys.argv[1] + ".tbl"

if int(sys.argv[2]) == 3:
    mpc_file = "input_data/" + sys.argv[1] + ".txt"
    wise_file = "input_data/" + sys.argv[1] + "_3band.tbl"

if int(sys.argv[2]) == 4:
    mpc_file = "input_data/" + sys.argv[1] + ".txt"
    wise_file = "input_data/" + sys.argv[1] + "_cryo.tbl"

source_ids = generate_source_ids_list(sys.argv[1], int(sys.argv[2]))
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