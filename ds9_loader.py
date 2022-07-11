"""
This program allows for certain FITS files to be loaded into DS9
"""
import sys
from mpc_wise_functions import comparer, generate_ra_dec
from ds9_interface_functions import make_region, load_files

load_file = sys.argv[1]
mpc_code = sys.argv[2]
bands = int(sys.argv[3])

load_files(load_file, mpc_code, bands)