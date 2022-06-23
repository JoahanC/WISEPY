#reference code from ds9_reader
from mpc_wise_functions import *
import os
import re


def make_region(file, w_band, source_ids):
    sid = file[:9:]
    with open(f"regions/{sid}_{w_band}.reg", 'w') as file:
        file.write("# Region file format: DS9 version 4.1\n")
        file.write('global color=green dashlist=8 3 width=1 ' 
                + 'font="helvetica 10 normal roman" select=1 ' 
                + 'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 ' 
                + 'include=1 source=1\n')
        file.write("fk5\n")
        file.write(f'circle({source_ids[sid][0]},{source_ids[sid][1]},15.000")')


mpc_file = "161989.txt"
wise_file = "table_irsa_catalog_search_results.tbl"
new_epochs = comparer(mpc_file, wise_file, True)

source_ids = []
for epoch in new_epochs:
    source_ids.append(new_epochs[epoch][0][:9])

files = os.listdir("WISE_Files/")
to_run = []
for file in files:
    if file[:9] in source_ids:
        to_run.append(file)

# Generating list of sorted 5 numbers keys
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
#print(key_sort)

# Sorting all files by first five numebers
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

# SORT BY LETTER PIECEWISE
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

#for file in w_sorted:
#    print(file)

# Renaming files to correct directory
renamed_sorted_run = []
for file in w_sorted:
    renamed_sorted_run.append("WISE_Files/" + file)

run_string = "ds9 -scale log -tile "
for file in renamed_sorted_run[:50]:
    run_string += file + ' -regions '
    reg_string = "regions/" + file[11:20] + '_' + file[21:23] + ".reg"
    run_string += reg_string + ' '
os.popen(run_string)