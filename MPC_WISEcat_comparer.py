#!/usr/bin/env python
# coding: utf-8
import csv 
import math
from astropy.time import Time
import pandas as pd
import numpy as np

# Utility functions


def decimal_day_converter(dec_day):
    """
    Turns a decimal day interpretation into a UTC format.
    Parameters: day
    Returns: String in the format: 'HH:MM:SS'
    """
    hour_remainder = float(dec_day) * 24
    hours = math.floor(hour_remainder)
    hours_str = str(hours).zfill(2)
    
    min_remainder = (hour_remainder - hours) * 60
    mins = math.floor(min_remainder)
    mins_str = str(mins).zfill(2)
    
    sec_remainder = (min_remainder - mins) * 60
    secs = sec_remainder
    secs_str = str(secs).zfill(2)
    
    
    return 'T' + hours_str + ':' + mins_str + ':' + secs_str


def MPC_parser(mpc_file):
    """
    Parses through the MPC datafile for a given object.
    Arguments: mpc_file (str) --  The file containing the relevant data
    Returns: (Np.Array) -- An array of date values for individual observations. 
    """

    read_file = pd.read_csv(mpc_file)
    read_file.to_csv ('mpc_file.csv', index=None)
    csv_url = 'mpc_file.csv'
    dates = []
    ids = []
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
            dec_day = day[2:]
            time_string = decimal_day_converter(dec_day)
            day = day[:2]
            date = year+'-'+month+'-'+day+time_string
            dates.append(date[:23])
    dates = np.array(dates)
    return dates


from pandas import *
def WISE_parser(wise_file):
    """
    Parses through the WISE datafile for a given object.
    Arguments: wise_file (str) -- The file containing the relevant data
    Returns: (Np.Array) -- An array of data values for individual observations.
    """
 
    data = read_csv(wise_file) 
    #setting up chart of data
    mjd = data['mjd']
    mjd_list = []
    for i in range(134):
        mjd_list.append(str(mjd[i]))
    #print(mjd_list)
    #Converting MJD to UTC
    times = mjd_list
    t = Time(times, format='mjd', scale='utc')

    #t.scale = 'utc'
    modified_time = t.isot
    #print(modified_time)

    #Adding UTC times to a dictionary
    Obs_Dates_dict = {}
    #adding mjd to dictionary
    Obs_Dates_dict["23606 1996 AS1"] = modified_time
    return modified_time


def base_round(x, base=5):
    return base * round(x/base)


# Comparison function between MPC and WISE catalog


def comparison(mpc_file, wise_file):
    """
    Compares the observational instances between the MPC and WISE dataset files for a given object. ds
    Arguments: mpc_file: txt file which contains MPC data
                wise_file: file which contains WISE
    """
    
    mpc_data = MPC_parser(mpc_file)
    wise_data = WISE_parser(wise_file)
    print("Raw Data \n")
    print("Epochs observed in the MPC database:", len(mpc_data))
    print("Epochs observed in the WISE database:", len(wise_data))
    
    # Year trim
    wise_years = []
    mpc_years = []
    trimmed_mpc_data = []
    for datum in wise_data:
        year = int(datum[:4])
        if year not in wise_years:
            wise_years.append(year)
    for datum in mpc_data:
        year = int(datum[:4])
        if year in wise_years:
            trimmed_mpc_data.append(datum)
        if year not in mpc_years:
            mpc_years.append(year)
    trimmed_mpc_data = np.array(trimmed_mpc_data)
    
    wise_years_str = ""
    for year in wise_years:
        wise_years_str += str(year) + ", "
    print("Years in which WISE data was collected: " + wise_years_str[:-2])
    mpc_years_str = ""
    for year in mpc_years:
        mpc_years_str += str(year) + ", "
    print("Years in which MPC data was collected: " + mpc_years_str[:-2])
    
    print("Epochs observed in trimmed MPC database:", len(trimmed_mpc_data))
    
    # Comparison sort
    rounded_mpc_data = []
    rounded_wise_data = []
    for datum in trimmed_mpc_data:
        pre_seconds = datum[:17]
        seconds = str(base_round(float(datum[17:]))).zfill(2)
        rounded_mpc_data.append(pre_seconds + seconds)
        
    for datum in wise_data:
        pre_seconds = datum[:17]
        seconds = str(base_round(float(datum[17:]))).zfill(2)
        rounded_wise_data.append(pre_seconds + seconds)
    
    new_epochs = np.array([epoch for epoch in rounded_wise_data if epoch not in rounded_mpc_data])
    print("New epochs observed in WISE catalog:", len(new_epochs))
    return rounded_mpc_data, rounded_wise_data, new_epochs


# Testing values 
mpc_data = '161989.txt'
wise_data = 'table_irsa_catalog_search_results-2.csv'
mpc_data, wise_data, new_epochs = comparison(mpc_data, wise_data)
#print(mpc_data)
#print(wise_data)
