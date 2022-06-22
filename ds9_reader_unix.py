from numpy import source
import mpc_wise_functions as mpc
import os


mpc_file = "161989.txt"
wise_file = "table_irsa_catalog_search_results.tbl"
new_epochs = mpc.comparer(mpc_file, wise_file, False)

source_ids = []
for epoch in new_epochs:
    source_ids.append(new_epochs[epoch][0][:9])


files = os.listdir("WISE_FILES/")
to_run = []
for file in files:
    if file[:9] in source_ids:
        to_run.append("WISE_FILES/" + file)


run_string = "ds9 -tile "
for file in to_run:
    run_string += file + ' '

os.popen(run_string)