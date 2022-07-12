"""
This program displays all available FITS files for viewing for a given
MPC designated object and a corresponding bandset.
"""
import os
import sys
from mpc_wise_functions import generate_source_ids_list
from ds9_interface_functions import generate_script


mpc_code = sys.argv[1]
bands = int(sys.argv[2])

uniques = len(generate_source_ids_list(mpc_code, bands))

view_all = input(f"There are {uniques} unique epochs detected.\n"
                 + "Would you like to look at all of them? [Y]/[N]\n")

if view_all.capitalize() == 'Y':

    script = generate_script(mpc_code, bands, 0, uniques * int(sys.argv[2]) + 1)
    
    os.popen(script)

else:

    lower_bound = input(f"Which image out of 1-{uniques * bands} would "
                        + "you want to start at\n")
    upper_bound = input(f"Which image out of 1-{uniques * bands} would "
                        + "you want to end at\n")
    script = generate_script(mpc_code, bands, int(lower_bound) - 1, 
                            int(upper_bound))
    os.popen(script)