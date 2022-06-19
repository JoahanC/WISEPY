import csv 
import math
from astropy.time import Time, TimeDelta
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
    Returns: (dict) -- the keys correspond to the utc dates for each image,
    a list is returned for each key where [0] corresponds to source_id and [1]
    corresponds to julian dates.
    """

    data_object = Table.read(wise_file, format='ipac')
    dates = list(data_object['mjd'])
    source_ids = list(data_object['source_id'])

    time_object = Time(dates, format='mjd', scale='utc')
    julian_dates = time_object.jd
    utc_dates = time_object.isot

    images = {}
    for idx, date in enumerate(utc_dates):
        images[str(date)] = [source_ids[idx], julian_dates[idx]]
    return images

def n_round(x, n=5):
    """
    Rounds a float to the nearest n integer.
    Arguments: x (float) -- the number being rounded.
               n (int) -- the number being rounded to.
    Returns: (int) -- the rounded integer.
    """
    return n * round(x/n)


def data_comparer(mpc_file, wise_file):
    """
    Compares the observational instances between the MPC and WISE dataset files for a given object. ds
    Arguments: mpc_file (str) -- txt file which contains MPC data
               wise_file (str) -- file which contains WISE
    """
    
    mpc_observation_data = MPC_parser(mpc_file)
    wise_image_data = WISE_parser(wise_file)
    print("-------------Database Statistics---------------- \n")
    print("Epochs observed in the MPC database:", len(mpc_observation_data))
    print("Epochs observed in the WISE database:", len(wise_image_data))
    
    # Year Counting
    wise_years = []
    mpc_years = []
    for datum in wise_image_data:
        print(datum)
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
    print("Years in which WISE data was collected: " + wise_years_str[:-2])
    mpc_years_str = ""
    for year in mpc_years:
        mpc_years_str += str(year) + ", "
    print("Years in which MPC data was collected: " + mpc_years_str[:-2])
    
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

    print(f"New epochs detected in WISE catalog: {len(new_epochs)}")
    return new_epochs