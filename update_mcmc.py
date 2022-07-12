"""
This file implements a script for generating input tables compatible with
tbl2MCMCin.py to be eventually used in Ned's MCMC model.
"""
import os
from ds9_interface_functions import generate_new_table, generate_full_table
from mpc_wise_functions import generate_flux_snr_plots
from epoch_JD_calendar import mjd_plot


files = os.listdir("loader_data")
epoch_files = {}
for file in files:
    if file.endswith("band.txt"):
        mpc_code = file[:len(file) - 10]
        if mpc_code in epoch_files.keys():
            epoch_files[mpc_code].append(file[:len(file)-8][-1:])
        else:
            epoch_files[mpc_code] = [file[:len(file)-8][-1:]]

band_lookup = {2: ".tbl", 3: "_3band.tbl", 4: "_cryo.tbl"}

for neo in epoch_files:
    for band in epoch_files[neo]:
        wise_file = f"input_data/{neo}{band_lookup[band]}"
        epoch_file = f"{neo}_{band}band.txt"
        print(f"Writing new input for {band} band data of {neo}.")
        generate_new_table(epoch_file, neo, band)
        print(f"Generating MJD plot for {band} band data of {neo}.")
        mjd_plot(neo, band)
        print(f"Writing MCMC input for {band} band data of {neo}.")
        generate_full_table(band, neo)
        print(f"Generating comparison plots for {band} band data of {neo}.")
        generate_flux_snr_plots(neo, band)