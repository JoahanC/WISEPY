import csv 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
from astropy.table import Table
from helpers import *


def MPC_parser(mpc_code):
    """
    Parses a txt file containing MPC data that returns utc and jd time values.
    
    Arguments: mpc_code (str) -- The MPC designated code for the object
    
    Returns: (dict) -- the keys correspond to the utc date for each observation,
    a list is returned for each key where [0] corresponds to the observation code
    and [1] corresponds to the julian date.
    """

    read_file = pd.read_csv(return_input_files(mpc_code)[0])
    read_file.to_csv('generated_data/mpc_file.csv', index=None)
    csv_url = 'generated_data/mpc_file.csv'
    dates = []
    ids = []
    obs_ids = []

    with open(csv_url) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            row_string = row[0]
            id_string = row_string[:15] 
            ids.append(id_string)
            row_string = row_string[15:]
            year = row_string[:4]
            month = row_string[5:7]
            row_string = row_string[8:]
            elements = row_string.split()
            day = elements[0]
            dec_day = day[2:9]
            time_string = decimal_day_converter(dec_day)
            day = day[:2]
            date = year+'-'+month+'-'+day+time_string
            dates.append(date[:23])
            observation_id = row_string[47:]
            obs_ids.append(observation_id)
    utc_dates = np.array(dates)
    time_object = Time(utc_dates, format='isot')
    julian_dates = time_object.jd
    observations = {}
    for idx, date in enumerate(utc_dates):
        observations[date] = [obs_ids[idx], julian_dates[idx]]
    return observations


def pull_fluxes(data_object, bands):
    """
    Pulls all flux values and flux sigmas for a given object.

    Arguments: data_object (Table) -- the data for a given object.
               bands (int) -- the bands included in the object data
    
    Returns: (tuple) -- All relevant flux values for the inputted data
    """
    w1_flux = list(data_object['w1flux'])
    w1_flux_sigma = list(data_object['w1sigflux'])
    w2_flux = list(data_object['w2flux'])
    w2_flux_sigma = list(data_object['w2sigflux'])
    if bands == 2:
        return (w1_flux, w1_flux_sigma,
               w2_flux, w2_flux_sigma)
    if bands == 3:
        w3_flux = list(data_object['w3flux'])
        w3_flux_sigma = list(data_object['w3sigflux'])
        return (w1_flux, w1_flux_sigma,
               w2_flux, w2_flux_sigma,
               w3_flux, w3_flux_sigma)
    if bands == 4:
        w3_flux = list(data_object['w3flux'])
        w3_flux_sigma = list(data_object['w3sigflux'])
        w4_flux = list(data_object['w4flux'])
        w4_flux_sigma = list(data_object['w4sigflux'])
        return (w1_flux, w1_flux_sigma,
               w2_flux, w2_flux_sigma,
               w3_flux, w3_flux_sigma,
               w4_flux, w4_flux_sigma)


def WISE_parser(mpc_code, bands=2):
    """
    Parses a tbl file containing WISE image data that returns source ids, utc,
    and jd time values.
    
    Arguments: mpc_code (str) -- The MPC designated code for the object
               bands (int) -- The targeted band set
    
    Returns: (dict) -- the keys correspond to the utc dates for each image,
    a list is returned for each key where [0] corresponds to source_id, [1]
    corresponds to julian dates, [2] corresponds to the ra, and [3] corresponds
    to nthe dec
    """

    data_object = Table.read(return_input_files(mpc_code, bands)[1], format='ipac')
    dates = list(data_object['mjd'])
    source_ids = list(data_object['source_id'])
    ra = list(data_object['ra'])
    dec = list(data_object['dec'])
    fluxes = pull_fluxes(data_object, bands)
    
    time_object = Time(dates, format='mjd', scale='utc')
    julian_dates = time_object.jd
    utc_dates = time_object.isot

    images = {}
    for idx, date in enumerate(utc_dates):
        if bands == 2:
            images[str(date)] = [source_ids[idx], julian_dates[idx], 
                                ra[idx], dec[idx], 
                                fluxes[0][idx], fluxes[1][idx], 
                                fluxes[2][idx], fluxes[3][idx]]
        if bands == 3:
            images[str(date)] = [source_ids[idx], julian_dates[idx], 
                                ra[idx], dec[idx], 
                                fluxes[0][idx], fluxes[1][idx], 
                                fluxes[2][idx], fluxes[3][idx],
                                fluxes[4][idx], fluxes[5][idx]]
        if bands == 4:
            images[str(date)] = [source_ids[idx], julian_dates[idx], 
                                ra[idx], dec[idx], 
                                fluxes[0][idx], fluxes[1][idx], 
                                fluxes[2][idx], fluxes[3][idx],
                                fluxes[4][idx], fluxes[5][idx],
                                fluxes[6][idx], fluxes[7][idx]]
    return images


def comparer(mpc_code, bands=2, show_stats=False):
    """
    Compares the observational instances between the MPC and WISE dataset 
    files for a given object. 
    
    Arguments: mpc_file (str) -- txt file which contains MPC data
               wise_file (str) -- file which contains WISE
    
    Returns: A dictionary containing all unique instances sorted by utc
    time with assorted information.
    """
    
    mpc_observation_data = MPC_parser(mpc_code)
    wise_image_data = WISE_parser(mpc_code, bands)
    
    # Year Counting
    wise_years = []
    mpc_years = []
    for datum in wise_image_data:
        year = int(datum[:4])
        if year not in wise_years:
            wise_years.append(year)
    for datum in mpc_observation_data:
        year = int(datum[:4])
        if year not in mpc_years:
            mpc_years.append(year)
    
    wise_years_str = ""
    for year in wise_years:
        wise_years_str += str(year) + ", "
    
    mpc_years_str = ""
    for year in mpc_years:
        mpc_years_str += str(year) + ", "
    
    # Acquire utc dates
    wise_utc = list(wise_image_data.keys())
    mpc_utc = list(mpc_observation_data.keys())

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

    new_epochs = {}
    for epoch in wise_intervals:
        recorded = False
        for observation in mpc_intervals:
            if wise_intervals[epoch] >= mpc_intervals[observation][0] and wise_intervals[epoch] <= mpc_intervals[observation][1]:
                recorded = True
        if not recorded:
            new_epochs[epoch] = wise_image_data[epoch]

    if show_stats:
        print("Epochs observed in the MPC database:", len(mpc_observation_data))
        print("Epochs observed in the WISE database:", len(wise_image_data))
        print("Years in which MPC data was collected: " + mpc_years_str[:-2])
        print("Years in which WISE data was collected: " + wise_years_str[:-2])
        print(f"New epochs detected in WISE catalog: {len(new_epochs)}")
        print('-' * 42)
        print(' ' * 7 + 'Time (utc)' + ' ' * 7 + '|' + ' ' * 4 + 'Source ID' + ' ' * 3)
        print('-' * 42)
        for epoch in new_epochs:
            print(f"{epoch}", '|', f"{new_epochs[epoch][0]}" )

    return new_epochs


def generate_source_ids_list(mpc_code, bands):
    """
    Generates a list of unique source_ids for a given object

    Arguments: mpc_code (str) -- The MPC designated code for the object
               bands (int) -- The targeted band set

    Returns: (list) -- A list of all unique source ids
    """
    new_epochs = comparer(mpc_code, bands, False)
    source_ids = []
    for epoch in new_epochs:
        source_ids.append(new_epochs[epoch][0][:9])
    return source_ids


def generate_ra_dec(mpc_code, bands):
    """
    Generates a dictionary of unique epochs sorted by source id and 
    containing the RA and DEC for the detected object in the corresponding
    FITS file.

    Arguments: mpc_code (str) -- The MPC designated code for the object
               bands (int) -- The targeted band set

    Returns: (dict) -- A dictionary of source_ids with the ids being the keys,
             the ra values[0], and the dec values[1]
    """

    data_object = Table.read(return_input_files(mpc_code, bands)[1], format='ipac')
    source_ids_2 = list(data_object['source_id'])
    ra = list(data_object['ra'])
    dec = list(data_object['dec'])
    info = {}
    for index, source_id in enumerate(source_ids_2):
        info[source_id[:9]] = [ra[index], dec[index]]

    return info


def generate_flux_snr_plots(mpc_code, band):
    """
    Generates a series of flux and SNR plots with associated error bars
    for a given MPC object and its targeted bandset.

    Arguments: mpc_code (str) -- The MPC designated code for the object
               bands (int) -- The targeted band set

    Returns: (None) -- A series of plots in /flux_plot and /snr_plot
    """
    new_data = Table.read(f"new_inputs/{mpc_code}_{band}bands.tbl", format='ipac')
    mjd_new = list(new_data['mjd'])
    w1_new = list(new_data['w1flux'])
    w2_new = list(new_data['w2flux'])
    w1_error_new = list(new_data["w1sigflux"])
    w2_error_new = list(new_data["w2sigflux"])
    snr1_new = list(new_data["w1snr"])
    snr2_new = list(new_data["w2snr"])

    all_data = Table.read(f"mcmc_inputs/{mpc_code}_{band}bands.tbl", format='ipac')
    mjd_all = list(all_data['mjd'])
    w1_all = list(all_data['w1flux'])
    w2_all = list(all_data['w2flux'])
    w1_error_all = list(all_data["w1sigflux"])
    w2_error_all = list(all_data["w2sigflux"])
    snr1_all = list(all_data["w1snr"])
    snr2_all = list(all_data["w2snr"])

    for i in range(len(mjd_new)):
        mjd_new[i] = float(mjd_new[i])
    for i in range(len(mjd_all)):
        mjd_all[i] = float(mjd_all[i])

    plt.errorbar(mjd_all, w1_all, label='Previous W1', marker=".", color="black", yerr=w1_error_all, fmt='o')
    plt.errorbar(mjd_new, w1_new, label='New W1', marker=".", color="red", yerr=w1_error_new, fmt='o')
    plt.title("Recorded flux epochs")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("Flux")
    plt.legend()
    plt.savefig(f"flux_plot/{mpc_code}_W1.png")
    plt.clf()

    plt.errorbar(mjd_all, w2_all, label='Previous W2', marker=".", color="black", yerr=w2_error_all, fmt='o')
    plt.errorbar(mjd_new, w2_new, label='New W2', marker=".", color="red", yerr=w2_error_new, fmt='o')
    plt.title(f"All Recorded Fluxes: {mpc_code}")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("Flux")
    plt.legend()
    plt.savefig(f"flux_plot/{mpc_code}_W2.png")
    plt.clf()

    plt.errorbar(mjd_new, w1_new, label='W1', marker=".", color="black", yerr=w1_error_new, fmt='o')
    plt.title(f"New Recorded Fluxes: {mpc_code}")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("Flux")
    plt.legend()
    plt.savefig(f"flux_plot/new_{mpc_code}_W1.png")
    plt.clf()

    plt.errorbar(mjd_new, w2_new, label='W2', marker=".", color="black", yerr=w2_error_new, fmt='o')
    plt.title(f"New Recorded Fluxes: {mpc_code}")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("Flux")
    plt.legend()
    plt.savefig(f"flux_plot/new_{mpc_code}_W2.png")
    plt.clf()

    plt.scatter(mjd_all, snr1_all, label='Previous W1', marker=".", color="black")
    plt.scatter(mjd_new, snr1_new, label='New W1', marker=".", color="red")
    plt.title(f"All Recorded SNRs: {mpc_code}")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("SNR")
    plt.legend()
    plt.savefig(f"snr_plot/{mpc_code}_W1.png")
    plt.clf()

    plt.scatter(mjd_new, snr1_new, label="W1", marker=".", color="red")
    plt.title(f"New Recorded SNRs: {mpc_code}")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("Flux")
    plt.legend()
    plt.savefig(f"snr_plot/new_{mpc_code}_W1.png")
    plt.clf()

    plt.scatter(mjd_all, snr2_all, label='Previous W2', marker=".", color="black")
    plt.scatter(mjd_new, snr2_new, label='New W2', marker=".", color="red")
    plt.title(f"All Recorded SNRs: {mpc_code}")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("SNR")
    plt.legend()
    plt.savefig(f"snr_plot/{mpc_code}_W2.png")
    plt.clf()

    plt.scatter(mjd_new, snr2_new, label="W2", marker=".", color="red")
    plt.title(f"New Recorded SNRs: {mpc_code}")
    plt.xlabel("Modified Julian Days")
    plt.ylabel("Flux")
    plt.legend()
    plt.savefig(f"snr_plot/new_{mpc_code}_W2.png")
    plt.clf()


