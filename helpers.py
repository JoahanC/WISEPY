"""
This file contains many computational utility functions used 
throughout the package.
"""
import math
import warnings
import matplotlib.pyplot as plt
import astropy
import datetime
import numpy as np
from matplotlib.ticker import MaxNLocator


_mpc_hex = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def return_input_files(packed_name, band=2):
    """
    Returns the valid MPC and WISE file for a given mpc object.

    Parameters
    ----------
    mpc_code : str
        The unpacked MPC designated code for the object.
    
    band : int
        The band set being queried.

    Returns
    -------
    The MPC file [0] and WISE file [1] associated with the object.
    """

    mpc_file = f"database_files/mpc/{packed_name}.txt"
    wise_file = f"database_files/wise/{packed_name}_{band}band.tbl"
    return mpc_file, wise_file


def decimal_day_converter(dec_day):
    """
    Turns a decimal day interpretation into a UTC format.
    
    Arguments
    ---------
    dec_day : str
        A representation of a fractional day; ie. 0.5, 0.90, 0.03, 0.11214, etc. 
    
    Returns
    -------
    A day interpretation in the format: 'HH:MM:SS'
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

    Parameters
    ----------
    x : float
        The number being rounded.
               
    n : int
        The number being rounded to.
    """
    return n * round(x/n)


def unpack_MPC_name(compact):
    if len(compact) == 5:
        if compact[4] == "P":
            # Periodic comet
            outn = "{:s}P".format(str(int(compact[0:4])))
        elif compact[4] == "S":
            # Natural Satellite
            outn = compact
        elif compact[0] == "~":
            # Numbered object at or above 620000
            dig1 = _mpc_hex.index(compact[1]) * (62**3)
            dig2 = _mpc_hex.index(compact[2]) * (62**2)
            dig3 = _mpc_hex.index(compact[3]) * 62
            dig4 = _mpc_hex.index(compact[4])
            num = 620000 + dig1 + dig2 + dig3 + dig4
            outn = "({:s})".format(str(num))
        else:
            # Numbered object
            outn = "({:s})".format(
                str(_mpc_hex.index(compact[0]) * 10000 + int(compact[1:]))
            )
    elif len(compact) == 7:
        # PLS object
        if compact[0:3] in ["PLS", "T1S", "T2S", "T3S"]:
            outn = "{:s} {:s}-{:s}".format(compact[3:], compact[0], compact[1])
        else:
            # Unnumbered asteroid
            year = int(_mpc_hex.index(compact[0]) * 100) + int(compact[1:3])
            if compact[4:6] == "00":
                prov = "{:1s}{:1s}".format(compact[3], compact[6])
            else:
                prov = "{:1s}{:1s}{:d}".format(
                    compact[3],
                    compact[6],
                    int(_mpc_hex.index(compact[4])) * 10 + int(compact[5]),
                )
            outn = "{:s} {:s}".format(str(year), prov)
    elif len(compact) == 8:
        # parabolic comet
        year = _mpc_hex.index(compact[1]) * 100 + int(compact[2:4])
        if compact[7] not in _mpc_hex[0:10]:
            prov = "{:s}{:s}{:s}".format(
                compact[4],
                compact[7],
                str(_mpc_hex.index(compact[5]) * 10 + int(compact[6])),
            )
            outn = "{:s}/{:s} {:s}".format(compact[0], str(year), prov)
        else:
            outn = "{:s}/{:s} {:s}{:s}".format(
                compact[0], str(year), compact[4], str(int(compact[5:7]))
            )
    else:
        raise NotImplementedError(
            "This designation could not be unpacked {:s}".format(compact)
        )
    return outn


def pack_MPC_name(name):
    """
    Converts an unpacked MPC designation to a packed MPC designation.

    Parameters
    ----------

    name : int
        An integer representation of an MPC name.

    Returns
    -------
    A string representation of a packed MPC name.
    """
    if name < 100000:
        return "%05d" % name

    elif name < 620000:
        nn = int(name / 10000.0)
        c = _mpc_hex[nn]
        return "{:1s}{:04d}".format(c, name - nn * 10000)

    else:
        # For numbers larger than 620,000 the MPC has defined a new packing
        # scheme. Code by J. Masiero
        nn = name - 620000
        dig4 = int(nn % 62)
        hold3 = (nn - dig4) / 62.0
        dig3 = int(hold3 % 62)
        hold2 = (hold3 - dig3) / 62.0
        dig2 = int(hold2 % 62)
        hold1 = (hold2 - dig2) / 62.0
        dig1 = int(hold1 % 62)
        return "~{:1s}{:1s}{:1s}{:1s}".format(
            _mpc_hex[dig1], _mpc_hex[dig2], _mpc_hex[dig3], _mpc_hex[dig4])


def template_new_plot(packed_name, mjd, flux_values, type, band, group):
    warnings.filterwarnings('error', module='astropy._erfa')
    fig, ax = plt.subplots()
    ax.errorbar(mjd, flux_values[0], yerr=flux_values[1], color="red", fmt=".", capsize=2)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.ticklabel_format(useOffset=False)
    ax.set_xlabel("Modified Julian Days")
    ax.set_ylabel(f"W{band} {type}")
    ax.set_title(f"{packed_name}", loc="left")
    fig.savefig(f"./plots/{type.lower()}_plots/{packed_name}/{group}_band/new_{type.lower()}_w{band}", dpi=1000)
    plt.close(fig)

def template_composite_plot(packed_name, new_mjd, all_mjd, new_flux_values, all_flux_values, type, band, group):
    warnings.filterwarnings('error', module='astropy._erfa')
    fig, ax = plt.subplots()
    ax.errorbar(all_mjd, all_flux_values[0], yerr=all_flux_values[1], color="black", fmt=".", capsize=2)
    ax.errorbar(new_mjd, new_flux_values[0], yerr=new_flux_values[1], color="red", fmt=".", capsize=2, ecolor="red")
    ax.ticklabel_format(useOffset=False)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Modified Julian Days")
    ax.set_ylabel(f"All W{band} {type}")
    ax.set_title(f"{packed_name}", loc="left")
    fig.savefig(f"./plots/{type.lower()}_plots/{packed_name}/{group}_band/all_{type.lower()}_w{band}", dpi=1000)
    plt.close(fig)


def terminal_table(mpc_code, band, good_epochs):
    """
    Generates an output table of source ids, frames indicies, and flux values,
    for a given set of epochs and a given band set.
    
    Arguments: mpc_code (str) -- the mpc designation for the asteroid
               band (int) -- the w band being targeted
               good_epochs (list) -- a list of valid source ids
    Returns: None (terminal printout)
    """

    band_file_map = {2: ("input_data/" + mpc_code + ".txt", 
                    "input_data/" + mpc_code + ".tbl", 0),
                    3: ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_3band.tbl", 22),
                    4: ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_cryo.tbl", 44)}

    # Table formatting below
    dash_string = '-' * 65 + '-' * band_file_map[band][2]
    test_string = "| Frame | Source Id "
    for i in range(1, band + 1):
        test_string += f"|  W{i} Flux | W{i} Sigma "
    test_string += '|'
    print(dash_string)
    print(test_string)
    print(dash_string)
    for idx, epoch in enumerate(good_epochs):
        epoch_string = '|' + f"{idx * 2 + 1}-{idx * 2 + 2}".rjust(6)
        epoch_string += f" | {epoch}"
        for i in range(int(band) * 2):
            epoch_string += ' | ' + f"{good_epochs[epoch][4 + i]}".rjust(8) 
        epoch_string += " | "
        print(epoch_string)
    print(dash_string)