from mpc_wise_functions import comparer, generate_source_ids_dict
from ds9_interface_functions import make_region, terminal_table
import os
import sys


def load_files(load_file, mpc_code, bands=2):
    file_stubs = []
    with open(f"loader_data/{load_file}", 'r') as file:
        for line in file:
            file_stubs.append(line.rstrip())

    wise_files = os.listdir(f"wise_images/{mpc_code}")
    sorted_run = []
    for stub in file_stubs:
        for file in wise_files:
            if file[:9] == stub:
                sorted_run.append(file)

    if bands == 4:
        idx = list(range(len(sorted_run)))[::4]
        w_sorted = []
        for i in idx:
            quartet = {int(sorted_run[i][11]) : sorted_run[i], 
               int(sorted_run[i+1][11]): sorted_run[i+1], 
               int(sorted_run[i+2][11]): sorted_run[i+2], 
               int(sorted_run[i+3][11]): sorted_run[i+3]}
            for i in range(1, 5):
                w_sorted.append(quartet[i])

    if bands == 3:
        idx = list(range(len(sorted_run)))[::3]
        w_sorted = []
        for i in idx:
            triplet = {int(sorted_run[i][11]) : sorted_run[i], 
               int(sorted_run[i+1][11]): sorted_run[i+1], 
               int(sorted_run[i+2][11]): sorted_run[i+2]}
            for i in range(1, 4):
                w_sorted.append(triplet[i])

    if bands == 2:
        idx = list(range(len(sorted_run)))[::2]
        w_sorted = []
        for i in idx:
            pair = {int(sorted_run[i][11]) : sorted_run[i], 
               int(sorted_run[i+1][11]): sorted_run[i+1]}
            for i in range(1, 3):
                w_sorted.append(pair[i])

    band_file_map = {2: ("input_data/" + mpc_code + ".txt", 
                    "input_data/" + mpc_code + ".tbl", 0),
                    3: ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_3band.tbl", 22),
                    4: ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_cryo.tbl", 44)}

    mpc_file, wise_file = band_file_map[bands][0], band_file_map[bands][1]

    new_epochs = comparer(mpc_code, bands, False)
    good_epochs = {}
    for epoch in new_epochs:
        sid = new_epochs[epoch][0][:9]
        for stub in file_stubs:
            if stub == sid:
                good_epochs[sid] = new_epochs[epoch]
                good_epochs[sid][0] = epoch    
  
    terminal_table(mpc_file, bands, good_epochs)
    
    region_list = []
    for file in w_sorted:
        region_list.append(f"wise_images/{mpc_code}/" + file)

    region_files = os.listdir("regions")
    for file in region_files:
        os.remove("regions/" + file)
    lookup_table = generate_source_ids_dict(mpc_file, wise_file, bands)
    for file in region_list:
        make_region(file, lookup_table, mpc_code)

    ds9_script = "ds9 -tile "
    for file in w_sorted:
        sid, w_band = file[:9], file[10:12]
        ds9_script += f"wise_images/{mpc_code}/" + file + ' -regions '
        ds9_script += f"regions/{sid}_{w_band}.reg "
    os.popen(ds9_script)


load_files(sys.argv[1], sys.argv[2], int(sys.argv[3]))