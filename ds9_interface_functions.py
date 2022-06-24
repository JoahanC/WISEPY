#reference code from ds9_reader
from mpc_wise_functions import *
import os
import re
import sys


def make_region(file, source_ids, mpc_code):
    """
    Makes a region file for a given .fits file in the regions folder.
    Arguments: file (str) -- the .fits file
               w_band (str) -- the corresponding w band for the file
               source_ids (dict) -- the source_id lookup table
    Returns: None (file in regions folder)
    """
    sid = file[13 + len(mpc_code):22 + len(mpc_code)]
    band = file[23 + len(mpc_code):25 + len(mpc_code)]
    with open(f"regions/{sid}_{band}.reg", 'w') as file:
        file.write("# Region file format: DS9 version 4.1\n")
        file.write('global color=green dashlist=8 3 width=1 ' 
                + 'font="helvetica 10 normal roman" select=1 ' 
                + 'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 ' 
                + 'include=1 source=1\n')
        file.write("fk5\n")
        file.write(f'circle({source_ids[sid][0]},{source_ids[sid][1]},15.000")')


def data_sort(source_ids, mpc_code):
    """
    The sorting algorithm for a given object's dataset.
    Arguments: source_ids (list) - a list of all unique epoch source ids.
    Returns: (list) - a sorted list of files corresponding to the source ids
             provided.
    """
    
    files = os.listdir("wise_images/" + mpc_code + '/')
    to_run = []
    for file in files:
        if file[:9] in source_ids:
            to_run.append(file)

    file_id = {}
    for file in to_run:
        if file == ".DS_Store" or file == ".DS_S":
            pass
        if int(file[:5]) in file_id.keys():
            file_id[int(file[:5])].append(file)
        else: 
            file_id[int(file[:5])] = [file]

    key_sort = list(file_id.keys())
    key_sort.sort()

    number_sorted = {}
    for key in key_sort:
        for file in to_run:
            code = int(file[0:5])
            if code == key and code in number_sorted.keys():
                number_sorted[code].append(file)
            elif code == key and code not in number_sorted.keys():
                number_sorted[code] = [file]
            else:
                pass

    for number in number_sorted:
        regex = r"[0-9]{5}[a]"
        matches = []
        for test in number_sorted[number]:
            if len(re.findall(regex, test)) == 0:
                pass
            else:
                matches.append(test)

        regex = r"[0-9]{5}[b]"
        for test in number_sorted[number]:
            if len(re.findall(regex, test)) == 0:
                pass
            else:
                matches.append(test)

        number_sorted[number] = matches

    # Putting sorted files into run list
    sorted_run = []
    for number in number_sorted:
        for file in number_sorted[number]:
            sorted_run.append(file)

    # W band sort
    idx = list(range(len(sorted_run)))[::2]

    w_sorted = []
    for i in idx:
        pair = (sorted_run[i], sorted_run[i + 1])
        if pair[0][10:12] == "w1":
            w_sorted.extend([pair[0], pair[1]])
        else:
            w_sorted.extend([pair[1], pair[0]])

    # Renaming files to correct directory
    renamed_sorted_run = []
    for file in w_sorted:
        renamed_sorted_run.append("wise_images/" + mpc_code + '/' + file)

    return renamed_sorted_run


def generate_script(source_ids, lower_bound, upper_bound, mpc_code):
    """
    Generates the ds9 script for running ds9_viewer.py.
    Arguments: source_ids (list) - A list of all relevant source ids to run
               lower_bound (int) - the initial index to load from the ids
               upperbound (int) - the ending index to load from the ids
    Returns: (str) - a set of ds9 commands to run in the terminal
    """
    sorted_files = data_sort(source_ids, mpc_code)

    region_files = os.listdir("regions")
    for file in region_files:
        os.remove("regions/" + file)

    mpc_file = "input_data/" + sys.argv[1] + ".txt"
    wise_file = "input_data/" + sys.argv[1] + ".tbl"

    lookup_table = generate_source_ids_dict(mpc_file, wise_file)

    for file in sorted_files:
        sid = file[13 + len(mpc_code):22 + len(mpc_code)]
        make_region(file, lookup_table, mpc_code)

    file_region = {}
    for file in sorted_files:
        sid = file[13 + len(mpc_code):22 + len(mpc_code)]
        band = file[23 + len(mpc_code):25 + len(mpc_code)]
        file_region[file] = f"regions/{sid}_{band}.reg"

    run_string = "ds9 -scale log -tile "
    for file in list(file_region.keys())[lower_bound:upper_bound]:
        run_string += file + ' -regions '
        reg_string = file_region[file]
        run_string += reg_string + ' '
    run_string += ' -zmax'
    return run_string