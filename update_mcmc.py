from ds9_interface_functions import write_epochs
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
        write_epochs(wise_file, epoch_file, neo, band)
