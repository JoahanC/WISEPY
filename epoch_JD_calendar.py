from visualizer_functions import visualizer
from mpc_wise_functions import *
import sys
from mpc_wise_functions import *
import matplotlib.pyplot as plt
import numpy as np
from astropy.utils.exceptions import AstropyUserWarning
import warnings
from astropy.time import Time


def mjd_plot(mpc_code, bands):
    
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    warnings.simplefilter('ignore', category=UserWarning)    

    #Setting up x-axis for histogram
    #band_lookup = {'2': "", '3': "_3band.tbl", '4': "_cryo.tbl"}

    mpc_file = "input_data/" + mpc_code + ".txt"
    #wise_file = "input_data/" + mpc_code + ".tbl"
    wise_file = f"new_inputs/{mpc_code}_{bands}bands.tbl"
    #print(wise_file)
    new_epochs = comparer(mpc_code, bands, False)

    JD_xaxis = []
    for epoch in new_epochs:
        jd_time = Time(new_epochs[epoch][1], format="jd")
        JD_xaxis.append(jd_time.mjd) #appending each JD for each epoch in the new epochs dictionary
   
    #Setting up y-axis for histogram

    
    JD_xaxis.sort()
    if len(JD_xaxis) != 0:
        plt.hist(JD_xaxis, bins=np.arange(min(JD_xaxis), max(JD_xaxis) + 1, 1), edgecolor="k", align='left', color='black')

        plt.title("New Epochs found per MJD")
        plt.xlabel("Julian Days")
        plt.ylabel("Unique Epochs")
        plt.savefig(f"mjd_plots/{mpc_code}_{bands}band.png")
        plt.close()