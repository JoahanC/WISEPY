from visualizer_functions import visualizer
from mpc_wise_functions import *
import sys
from mpc_wise_functions import *
import matplotlib.pyplot as plt
import numpy as np


def JD_plot(mpc_code, bands):
    #Setting up x-axis for histogram
    band_lookup = {'2': ".tbl", '3': "_3band.tbl", '4': "_cryo.tbl"}

    mpc_file = "input_data/" + mpc_code + ".txt"
    #wise_file = "input_data/" + mpc_code + ".tbl"
    wise_file = f"input_data/{mpc_code}{band_lookup[bands]}"

    new_epochs = comparer(mpc_file, wise_file, False)

    JD_xaxis = []
    for epoch in new_epochs:
        JD_xaxis.append(new_epochs[epoch][1]) #appending each JD for each epoch in the new epochs dictionary
   
    #Setting up y-axis for histogram

    
    JD_xaxis.sort()

    plt.hist(JD_xaxis, bins = 400, edgecolor="k", align='left', color='black')

    plt.title("New Epochs found per JD")
    plt.xlabel("Julian Days")
    plt.ylabel("Unique Epochs")
    plt.savefig(f"jd_plots/{mpc_code}_{bands}bands.png")
    #plt.show()
    plt.close()

#JD_plot(new_epochs, 'placeholder', mpc_code)