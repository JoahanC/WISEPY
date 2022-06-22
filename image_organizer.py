#reference code from ds9_reader
from unittest import skip
from mpc_wise_functions import *
import os


mpc_file = "161989.txt"
wise_file = "table_irsa_catalog_search_results.tbl"
new_epochs = comparer(mpc_file, wise_file, False)

source_ids = []
for epoch in new_epochs:
    source_ids.append(new_epochs[epoch][0][:9])

files = os.listdir("WISE_Files/")
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

# LETTER SORT
for code in file_id:
    unsorted = file_id[code]
    sorted = []
    for file in unsorted:
        if file[5:6] == 'a':
            sorted.append(file)
    for file in unsorted:
        if file[5:6] == 'b':
            sorted.append(file)
    file_id[code] = sorted

for vals in file_id:
    print(vals, file_id[vals]) 


flag_dict = {}
for file in files:
    flag_dict[file] = False

ordered_files = []
for file in files:
    if not flag_dict[file]:
        sid = file[:9]
        dupe_list = []
        for file in files:
            if file[:9] == sid:
                dupe_list.append(file)
        #print(len(dupe_list))
    




"""
image_list = []

for sid in source_ids:
    if sid != files:
        skip
    if sid == files:
        image_list.append(sid)
        if len(image_list) == 2:
            if sid[10:12] == 'w2':
                image_list.append(image_list.pop(sid))
            if sid[10:12] == 'w1':
                image_list.pop(sid)
                image_list.insert(0, sid)"""

