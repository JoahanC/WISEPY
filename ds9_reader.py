from numpy import source
import mpc_wise_functions as mpc
import os


mpc_file = "161989.txt"
wise_file = "table_irsa_catalog_search_results.tbl"
new_epochs = mpc.comparer(mpc_file, wise_file, False)

source_ids = {}
for epoch in new_epochs:
    source_ids[new_epochs[epoch][0][:9]] = [new_epochs[epoch][2], new_epochs[epoch][3]]

files = os.listdir("WISE_Files/")
to_run = []
for file in files:
    if file[:9] in source_ids:
        to_run.append("WISE_Files/" + file)

def make_region(file, w_band, source_ids):
    sid = file[:9:]
    with open(f"regions/{sid}_{w_band}.reg", 'w') as file:
        file.write("# Region file format: DS9 version 4.1\n")
        file.write('global color=green dashlist=8 3 width=1' 
                + 'font="helvetica 10 normal roman" select=1' 
                + 'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1' 
                + 'include=1 source=1\n')
        file.write("fk5\n")
        file.write(f'circle({source_ids[sid][0]},{source_ids[sid][1]},30.000")')

files2 = os.listdir("regions")
for f in files2:
    os.remove("regions/" + f)
counter = 0
for file in files:
    if file[:9] in source_ids:
        counter += 1
        make_region(file, file[10:12], source_ids)
print(counter)

"""
os.popen("ds9")
run_string = "open "
for file in renamed_sorted_run:
    run_string += file + ' '
os.popen(run_string)
os.popen("-single")"""