#reference code from ds9_reader
from mpc_wise_functions import *
import os
import re
import warnings
from astropy.time import Time, TimeDelta
from astropy.utils.exceptions import AstropyUserWarning
from terminal_helpers import terminal_table


def make_region(mpc_code, fits_file, source_ids):
    """
    Makes a region file for a given .fits file in the regions folder.
    Arguments: file (str) -- the .fits file
               w_band (str) -- the corresponding w band for the file
               source_ids (dict) -- the source_id lookup table
    Returns: None (file in regions folder)
    """
    sid = fits_file[13 + len(mpc_code):22 + len(mpc_code)]
    band = fits_file[23 + len(mpc_code):25 + len(mpc_code)]
    with open(f"regions/{sid}_{band}.reg", 'w') as file:
        file.write("# Region file format: DS9 version 4.1\n")
        file.write('global color=green dashlist=8 3 width=1 ' 
                + 'font="helvetica 10 normal roman" select=1 ' 
                + 'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 ' 
                + 'include=1 source=1\n')
        file.write("fk5\n")
        file.write(f'circle({source_ids[sid][0]},{source_ids[sid][1]},15.000")')


def data_sort(source_ids, mpc_code, bands=2):
    """
    The sorting algorithm for a given object's dataset.
    
    Arguments: source_ids (list) - a list of all unique epoch source ids.
    
    Returns: (list) - a sorted list of files corresponding to the source ids
             provided.
    """
    
    # Generate all relevant keys and perform initial sort
    files = os.listdir("wise_images/" + mpc_code + '/')
    wise_files = []
    for file in files:
        if file[:9] in source_ids:
            wise_files.append(file)

    file_id = {}
    for file in wise_files:
        if file == ".DS_Store" or file == ".DS_S":
            pass
        if int(file[:5]) in file_id.keys():
            file_id[int(file[:5])].append(file)
        else: 
            file_id[int(file[:5])] = [file]

    key_sort = list(file_id.keys())
    key_sort.sort()

    # Sky region sort
    number_sorted = {}
    for key in key_sort:
        for file in wise_files:
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
        regex = r"[0-9]{5}[c]"
        for test in number_sorted[number]:
            if len(re.findall(regex, test)) == 0:
                pass
            else:
                matches.append(test)
        regex = r"[0-9]{5}[r]"
        for test in number_sorted[number]:
            if len(re.findall(regex, test)) == 0:
                pass
            else:
                matches.append(test)
        regex = r"[0-9]{5}[s]"
        for test in number_sorted[number]:
            if len(re.findall(regex, test)) == 0:
                pass
            else:
                matches.append(test)

        number_sorted[number] = matches

    sorted_run = []
    for number in number_sorted:
        a_temp = {}
        b_temp = {}
        for elem in number_sorted[number]:
            if elem[5] == 'a':
                if elem[6:9] not in a_temp.keys():
                    a_temp[elem[6:9]] = [elem]
                else:
                    a_temp[elem[6:9]].append(elem)
            else:
                if elem[6:9] not in b_temp.keys():
                    b_temp[elem[6:9]] = [elem]
                else:
                    b_temp[elem[6:9]].append(elem)
        a_keys = list(a_temp.keys())
        b_keys = list(b_temp.keys())
        a_keys.sort()
        b_keys.sort()
        for key in a_keys:
            for elem in a_temp[key]:
                sorted_run.append(elem)
        for key in b_keys:
            for elem in b_temp[key]:
                sorted_run.append(elem)

    # Test W band sort

    idx = list(range(len(sorted_run)))[::bands]
    w_sorted = []
    interval = {}
    for i in idx:
        for j in range(bands):
            interval[int(sorted_run[i + j][11])] = sorted_run[i + j]
        for k in range(1, bands + 1):
            w_sorted.append(interval[k])

    renamed_sorted_run = []
    for file in w_sorted:
        renamed_sorted_run.append("wise_images/" + mpc_code + '/' + file)
    return renamed_sorted_run


def generate_script(mpc_code, bands, lower_bound, upper_bound):
    """
    Generates the ds9 script for running ds9_viewer.py.
    Arguments: source_ids (list) - A list of all relevant source ids to run
               lower_bound (int) - the initial index to load from the ids
               upperbound (int) - the ending index to load from the ids
    Returns: (str) - a set of ds9 commands to run in the terminal
    """

    source_ids = generate_source_ids_list(mpc_code, bands)
    sorted_files = data_sort(source_ids, mpc_code, bands)
    region_files = os.listdir("regions")
    for file in region_files:
        os.remove("regions/" + file)

    lookup_table = generate_ra_dec(mpc_code, bands)

    for file in sorted_files:
        sid = file[13 + len(mpc_code):22 + len(mpc_code)]
        make_region(mpc_code, file, lookup_table)

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


def generate_new_table(sid_file, mpc_code, band):
    """
    Generates an ipac format table from an existing WISE table for a set
    cluster of images from the IRSA WISE Image Service.
    Arguments: sid_file (str) -- the file name of the desired source ids
               mpc_code (str) -- the mpc code of the asteroid being observed
               band (str) -- the band being observed for the set of source ids
    Returns: None (returns new tbl files in mcmc_inputs)
    """
    
    warnings.simplefilter('ignore', category=AstropyUserWarning)

    # Reads in all WISE source ids for a given asteroid and band set
    band_lookup = {2: ".tbl", 3: "_3band.tbl", 4: "_cryo.tbl"}
    wise_file = f"input_data/{mpc_code}{band_lookup[band]}"
    data_object = Table.read(wise_file, format='ipac')
    wise_sids = list(data_object['source_id'])
    wise_sids = [sid[0:9] for sid in wise_sids]

    # Reads in unique source ids
    unique_sids = []
    with open(f"loader_data/{sid_file}", 'r') as file:
        for line in file:
            unique_sids.append(line.rstrip())
    # Generates a mask with unique source ids and writes tbl file
    mask = np.isin(wise_sids, unique_sids)
    t_new = data_object[mask]
    t_new.write(f"new_inputs/{mpc_code}_{band}bands.tbl", 
            format="ipac", overwrite=True)


def return_new_sids(band, mpc_code):
    """
    Generates an ipac format table from an existing WISE table for a set
    cluster of images from the IRSA WISE Image Service.
    Arguments: sid_file (str) -- the file name of the desired source ids
               mpc_code (str) -- the mpc code of the asteroid being observed
               band (str) -- the band being observed for the set of source ids
    Returns: None (returns new tbl files in mcmc_inputs)
    """
    
    warnings.simplefilter('ignore', category=AstropyUserWarning)

    # Reads in all WISE source ids for a given asteroid and band set
    new_sids = []
    try:
        wise_file = f"new_inputs/{mpc_code}_{band}bands.tbl"
        data_object = Table.read(wise_file, format='ipac')
        sids = list(data_object['source_id'])
        sids = [sid[0:9] for sid in sids]
        new_sids.extend(sids)
    except:
        print('')
    return new_sids


def generate_full_table(mpc_code, band):
    """
    Generates a TBL file containing all useful WISE observational instances
    which are useful for thermophysical modeling.

    Arguments: mpc_code (str) -- The MPC designated code for the object
               bands (int) -- The targeted band set

    Returns: (None) -- A new TBL file in /mcmc_inputs
    """
    warnings.simplefilter('ignore', category=AstropyUserWarning)

    # Reads in all known WISE observations
    band_lookup = {2: ".tbl", 3: "_3band.tbl", 4: "_cryo.tbl"}
    listy = []
    wise_file = f"input_data/{mpc_code}{band_lookup[band]}"
    try:
        listy.append(WISE_parser(mpc_code, int(band)))
    except:
        print('')
    
    wise_obs = {}
    for dict in listy:
        for date in dict:
            wise_obs[date] = dict[date]

    # Read in all known MPC observations
    mpc_file = "input_data/" + mpc_code + ".txt"
    mpc_obs = MPC_parser(mpc_code)
    reported_obs = {}
    for date in mpc_obs:
        if mpc_obs[date][0][7:] == "C51":
            reported_obs[date] = mpc_obs[date]

    # Acquire utc dates
    wise_utc = list(wise_obs.keys())
    mpc_utc = list(mpc_obs.keys())
    
    mpc_intervals = {}
    utc_delta = TimeDelta("11.0", format='sec')
    for datum in mpc_utc:
        base_time = Time(datum, format='isot')
        lower_bound = base_time - utc_delta
        upper_bound = base_time + utc_delta
        mpc_intervals[datum] = [lower_bound.jd, upper_bound.jd]

    wise_intervals = {}
    for datum in wise_utc:
        base_time = Time(datum, format='isot')
        wise_intervals[datum] = [base_time.jd]

    existing_epochs = {}
    for epoch in wise_intervals:
        recorded = False
        for observation in mpc_intervals:
            if wise_intervals[epoch] >= mpc_intervals[observation][0] and wise_intervals[epoch] <= mpc_intervals[observation][1]:
                recorded = True
        if recorded:
            existing_epochs[epoch] = wise_obs[epoch]

    existing_sids = []
    for epoch in existing_epochs:
        existing_sids.append(existing_epochs[epoch][0][:9])

    new_sids = return_new_sids(band, mpc_code)
    all_sids = existing_sids + new_sids

    band_lookup = {2: ".tbl", 3: "_3band.tbl", 4: "_cryo.tbl"}
    wise_file = f"input_data/{mpc_code}{band_lookup[band]}"
    data_object = Table.read(wise_file, format='ipac')
    wise_sids = list(data_object['source_id'])
    wise_sids = [sid[0:9] for sid in wise_sids]
    mask = np.isin(wise_sids, all_sids)
    t_new = data_object[mask]
    t_new.write(f"mcmc_inputs/irsa_{mpc_code}_{band}bands.tbl", 
            format="ipac", overwrite=True)


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
    lookup_table = generate_ra_dec(mpc_code, bands)
    for file in region_list:
        make_region(mpc_code, file, lookup_table)

    ds9_script = "ds9 -tile "
    for file in w_sorted:
        sid, w_band = file[:9], file[10:12]
        ds9_script += f"wise_images/{mpc_code}/" + file + ' -regions '
        ds9_script += f"regions/{sid}_{w_band}.reg "
    os.popen(ds9_script)