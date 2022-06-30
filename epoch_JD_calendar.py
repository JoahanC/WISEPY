from visualizer_functions import visualizer
from mpc_wise_functions import *
import sys
from mpc_wise_functions import *
import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
import numpy as np

#Obtaining epoch data
mpc_file = "input_data/" + sys.argv[1] + ".txt"
wise_file = "input_data/" + sys.argv[1] + ".tbl"
mpc_code = sys.argv[1]

new_epochs = comparer(mpc_file, wise_file, False)
#visualizer(mpc_file, wise_file, sys.argv[2])
#print(new_epochs)


def JD_plot(new_epochs, title, mpc_code):
    #Setting up x-axis for histogram

    JD_xaxis = []
    for epoch in new_epochs:
        JD_xaxis.append(new_epochs[epoch][1]) #appending each JD for each epoch in the new epochs dictionary
   
    #Setting up y-axis for histogram

    
    JD_xaxis.sort()
    for jd in JD_xaxis:
        print(jd)
    plt.hist(JD_xaxis, bins = 400, edgecolor="k", align='left', color='black')

    plt.title("New Epochs found per JD")
    plt.xlabel("Julian Days")
    plt.ylabel("Unique Epochs")
    plt.savefig(f"jd_plots/{mpc_code}.png")
    #plt.show()
    plt.close()

JD_plot(new_epochs, 'placeholder', mpc_code)