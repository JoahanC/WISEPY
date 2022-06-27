from mpc_wise_functions import comparer
import os
import sys


def load_files(load_file, mpc_code):
    file_stubs = []
    with open(f"loader_data/{load_file}", 'r') as file:
        for line in file:
            file_stubs.append(line.rstrip())

    wise_files = os.listdir(f"wise_images/{mpc_code}")
    run_files = []
    for stub in file_stubs:
        for file in wise_files:
            if file[:9] == stub:
                run_files.append(file)

    mpc_file = "input_data/" + mpc_code + ".txt"
    wise_file = "input_data/" + mpc_code + ".tbl"
    new_epochs = comparer(mpc_file, wise_file, False)
    good_epochs = {}
    for epoch in new_epochs:
        sid = new_epochs[epoch][0][:9]
        for stub in file_stubs:
            if stub == sid:
                good_epochs[sid] = new_epochs[epoch]
                good_epochs[sid][0] = epoch    
    print('-' * 57)
    print("| Source Id | W1 Flux  | W1 Sigma | W2 Flux  | W2 Sigma |")
    print('-' * 57)
    for epoch in good_epochs:
        print(f"| {epoch}" + ' | ' + f"{good_epochs[epoch][4]}".rjust(8) 
              + ' | ' + f"{good_epochs[epoch][5]}".rjust(8) 
              + ' | ' + f"{good_epochs[epoch][6]}".rjust(8)
              + ' | ' + f"{good_epochs[epoch][7]}".rjust(8) + ' | ')

    ds9_script = "ds9 -tile "
    for file in run_files:
        sid, w_band = file[:9], file[10:12]
        ds9_script += f"wise_images/{mpc_code}/" + file + ' -regions '
        ds9_script += f"regions/{sid}_{w_band}.reg "
    os.popen(ds9_script)


load_files(sys.argv[1], sys.argv[2])