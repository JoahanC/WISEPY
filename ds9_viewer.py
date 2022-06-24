from mpc_wise_functions import *
from ds9_interface_functions import *
import os
import re
import sys


mpc_file = "input_data/" + sys.argv[1] + ".txt"
wise_file = "input_data/" + sys.argv[1] + ".tbl"

source_ids = generate_source_ids(mpc_file, wise_file)

script = generate_script(source_ids)
os.popen(script)

