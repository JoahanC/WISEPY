import csv 
import math
from astropy.time import Time
from astropy.table import Table
import pandas as pd
import numpy as np
from numpy import source


def decimal_day_converter(dec_day):
    """
    Turns a decimal day interpretation into a UTC format.
    Parameters: dec_day (str) -- a representation of a fractional day; ie. 0.5,
    0.90, 0.03, 0.11214, etc. 
    Returns: (str) a day interpretation in the format: 'HH:MM:SS'
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
    Parses a txt file containing MPC data that returns utc and jd time values.
    
    Arguments: mpc_file (str) -- a path to the data file to read from.
    Returns: (dict) -- the keys correspond to the utc date for each observation,
    a list is returned for each key where [0] corresponds to the observation code
    and [1] corresponds to the julian date.
    """

    read_file = pd.read_csv(mpc_file)
    read_file.to_csv ('mpc_file.csv', index=None)
    csv_url = 'mpc_file.csv'
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
            dec_day = day[2:]
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


def WISE_parser(wise_file):
    """
    Parses a tbl file containing WISE image data that returns source ids, utc,
    and jd time values.
    
    Arguments: wise_file (str) -- a path to the data file to read from.
    Returns: (dict) -- the keys correspond to the source ids for each image,
    a list is returned for each key where [0] corresponds to jd and [1]
    corresponds to utc.
    """

    data_object = Table.read(wise_file, format='ipac')
    dates = list(data_object['mjd'])
    source_ids = list(data_object['source_id'])

    time_object = Time(dates, format='mjd', scale='utc')
    julian_dates = time_object.jd
    utc_dates = time_object.isot

    images = {}
    for idx, id in enumerate(source_ids):
        images[id] = [julian_dates[idx], utc_dates[idx]]

    return images

def n_round(x, n=5):
    """
    Rounds a float to the nearest n integer.
    Arguments: x (float) -- the number being rounded.
               n (int) -- the number being rounded to.
    Returns: (int) -- the rounded integer.
    """
    return n * round(x/n)


def comparison(mpc_file, wise_file):
    """
    Compares the observational instances between the MPC and WISE dataset files for a given object. ds
    Arguments: mpc_file: txt file which contains MPC data
                wise_file: file which contains WISE
    """
    
    mpc_data = MPC_parser(mpc_file)
    wise_data = WISE_parser(wise_file)
    print("-------------Database Statistics \n")
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
        #if year in wise_years:
        #    trimmed_mpc_data.append(datum)
        if year not in mpc_years:
            mpc_years.append(year)
    
    wise_years_str = ""
    for year in wise_years:
        wise_years_str += str(year) + ", "
    print("Years in which WISE data was collected: " + wise_years_str[:-2])
    mpc_years_str = ""
    for year in mpc_years:
        mpc_years_str += str(year) + ", "
    print("Years in which MPC data was collected: " + mpc_years_str[:-2])
    
    # Comparison sort
    rounded_mpc_data = []
    rounded_wise_data = []
    for datum in mpc_data:
        pre_seconds = datum[:17]
        seconds = str(n_round(float(datum[17:]))).zfill(2)
        rounded_mpc_data.append(pre_seconds + seconds)
        
    for datum in wise_data:
        pre_seconds = datum[:17]
        seconds = str(n_round(float(datum[17:]))).zfill(2)
        rounded_wise_data.append(pre_seconds + seconds)
    
    new_epochs = np.array([epoch for epoch in rounded_wise_data if epoch not in rounded_mpc_data])
    print("New epochs observed in WISE catalog:", len(new_epochs))
    return rounded_mpc_data, rounded_wise_data, new_epochs




