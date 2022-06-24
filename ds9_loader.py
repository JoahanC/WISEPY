import os
import sys


file_stubs = []
with open("loader_data/load.txt", 'r') as file:
    for line in file:
        file_stubs.append(line.rstrip())

wise_files = os.listdir(f"wise_images/{sys.argv[1]}")
run_files = []
for stub in file_stubs:
    for file in wise_files:
        if file[:9] == stub:
            run_files.append(file)

ds9_script = "ds9 -tile "
for file in run_files:
    sid, w_band = file[:9], file[10:12]
    ds9_script += f"wise_images/{sys.argv[1]}/" + file + ' -regions '
    ds9_script += f"regions/{sid}_{w_band}.reg "
os.popen(ds9_script)