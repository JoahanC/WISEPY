"""
This file contains many computational utility functions used 
throughout the package.
"""
import math


def return_input_files(mpc_code, bands=2):
    """
    Returns the valid MPC and WISE file for a given mpc object

    Arguments: mpc_code (str) -- the MPC designated code for an object
               bands (int) -- the number of bands associated with an 
               object

    Returns: (tuple) The MPC file [0] and WISE file [1] associated with
             the object
    """

    band_lookup_table = {2: ".tbl", 3: "_3band.tbl", 4 : "_cryo.tbl"}
    mpc_file = f"input_data/{mpc_code}.txt"
    wise_file = f"input_data/{mpc_code}{band_lookup_table[bands]}"
    return mpc_file, wise_file


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


def n_round(x, n=5):
    """
    Rounds a float to the nearest n integer.
    Arguments: x (float) -- the number being rounded.
               n (int) -- the number being rounded to.
    Returns: (int) -- the rounded integer.
    """
    return n * round(x/n)