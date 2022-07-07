"""
This file implements a script for generating input tables compatible with
tbl2MCMCin.py to be eventually used in Ned's MCMC model.
"""
from ds9_interface_functions import writeMCMC_table, generate_full_table
from mpc_wise_functions import flux_scatter
from epoch_JD_calendar import JD_plot
import warnings
import os


files = os.listdir("loader_data")
epoch_files = {}
for file in files:
    if file.endswith("band.txt"):
        mpc_code = file[:len(file) - 10]
        if mpc_code in epoch_files.keys():
            epoch_files[mpc_code].append(file[:len(file)-8][-1:])
        else:
            epoch_files[mpc_code] = [file[:len(file)-8][-1:]]

band_lookup = {'2': ".tbl", '3': "_3band.tbl", '4': "_cryo.tbl"}

for neo in epoch_files:
    for band in epoch_files[neo]:
        wise_file = f"input_data/{neo}{band_lookup[band]}"
        epoch_file = f"{neo}_{band}band.txt"
        writeMCMC_table(epoch_file, neo, band)
        JD_plot(neo, band)
        generate_full_table(band, neo)
        flux_scatter(neo, band)
        print(f"Writing MCMC input for {band} band data of {neo}.")
