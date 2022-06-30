from astropy.table import Table
from nbformat import write
from mpc_wise_functions import comparer
from ds9_interface_functions import data_sort
import os
import sys


def terminal_table(band_file_map, bands, good_epochs):

    dash_string = '-' * 65 + '-' * band_file_map[bands][2]
    test_string = "| Frame | Source Id "
    for i in range(1, int(bands) + 1):
        test_string += f"|  W{i} Flux | W{i} Sigma "
    test_string += '|'
    print(dash_string)
    print(test_string)
    print(dash_string)
    for idx, epoch in enumerate(good_epochs):
        epoch_string = '|' + f"{idx * 2 + 1}-{idx * 2 + 2}".rjust(6)
        epoch_string += f" | {epoch}"
        for i in range(int(bands) * 2):
            epoch_string += ' | ' + f"{good_epochs[epoch][4 + i]}".rjust(8) 
        epoch_string += " | "
        print(epoch_string)
    print(dash_string)

def load_epochs(load_file, mpc_code, bands):

    file_stubs = []
    with open(f"loader_data/{load_file}", 'r') as file:
        for line in file:
            file_stubs.append(line.rstrip())

    band_file_map = {'2': ("input_data/" + mpc_code + ".txt", 
                    "input_data/" + mpc_code + ".tbl", 0),
                    '3': ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_3band.tbl", 22),
                    '4': ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_cryo.tbl", 44)}

    mpc_file, wise_file = band_file_map[bands][0], band_file_map[bands][1]
    
    new_epochs = comparer(mpc_file, wise_file, False, int(bands))
    good_epochs = {}
    for epoch in new_epochs:
        sid = new_epochs[epoch][0][:9]
        for stub in file_stubs:
            if stub == sid:
                good_epochs[sid] = new_epochs[epoch]
                good_epochs[sid][0] = epoch    

    terminal_table(band_file_map, bands, good_epochs)
    return file_stubs, good_epochs



def write_table(load_file, band_file_map, bands, mpc_code):

    file_stubs = []
    with open(f"loader_data/{load_file}", 'r') as file:
        for line in file:
            file_stubs.append(line.rstrip())

    mpc_file, wise_file = band_file_map[bands][0], band_file_map[bands][1]
    new_epochs = comparer(mpc_file, wise_file, False, int(bands))

    good_epochs = {}
    for epoch in new_epochs:
        sid = new_epochs[epoch][0][:9]
        for stub in file_stubs:
            if stub == sid:
                good_epochs[sid] = new_epochs[epoch]
                good_epochs[sid][0] = epoch   

    sids = list(good_epochs.keys())
    for sid in sids:
        print(sid, good_epochs[sid])
    
    big_list = []
    for i in range(int(bands) * 2):
        temp_list = []
        for sid in sids:
            temp_list.append(good_epochs[sid][4+i])
        big_list.append(temp_list)
    

    idx_lookup = {-1: "source_id", 0: "w1flux", 1: "w1sigflux",
                  2: "w2flux", 3: "w2sigflux",
                  4: "w3flux", 5: "w3sigflux",
                  6: "w4flux", 7: "w4sigflux"}
    write_dict = {}
    write_dict[idx_lookup[-1]] = sids

    for i in range(int(bands) * 2):
        write_dict[idx_lookup[i]] = []
        for sid in sids:
            write_dict[idx_lookup[i]].append(good_epochs[sid][4+i])
    for key in write_dict:
        print(key, write_dict[key])

    for sid in good_epochs:
        print(sid, good_epochs[sid])

    first = list(write_dict.keys())
    last = list(write_dict.values())

    dat = Table(last, names=first)
    dat.write(f"mcmc_inputs/{mpc_code}_{bands}bands.tbl", 
          format="ipac", overwrite=True)

