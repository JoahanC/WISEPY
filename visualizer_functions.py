from mpc_wise_functions import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import datetime
import os


def month_plot(dates, title, title_year):
    """
    Produces a histogram of recorded epochs during a given year
    across different months.
    Arguments: dates (list) -- a list of utc format strings 
               title (str) -- the name of the dataset
               title_year (str) -- the year being plotted for. If
               a cumulative graph is made, provide the string:
               "all years"
    Returns: a histogram in the output folder
    """

    trimmed_dates = []
    for date in dates:
        year, month, day = date[:4], date[5:7], date[8:10]
        date_string = year + month + day
        trimmed_dates.append(date_string)

    dtm = lambda x: int(x[4:6])
    months = list(map(dtm, trimmed_dates))

    fig, ax = plt.subplots()
    bins = np.arange(1,14)
    ax.hist(months, bins = bins, edgecolor="k", align='left', color='black')
    ax.set_xticks(bins[:-1])
    ax.set_title(title)
    ax.set_xticklabels([datetime.date(1900,i,1).strftime('%b') for i in bins[:-1]])
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    if title == "Unique":
        plt.title(f"New epochs found during {title_year}")
        if title_year == "all years":
            title_year = "Cumulative"
        plt.xlabel("Month")
        plt.ylabel("Count")
        plt.savefig(f"output_plots/new_epochs_{title_year}")
        plt.close()
    else:
        plt.title(f"Epochs in the {title} dataset during {title_year}")
        if title_year == "all years":
            title_year = "Cumulative"
        plt.xlabel("Month")
        plt.ylabel("Count")
        plt.savefig(f"output_plots/{title}_{title_year}")
        plt.close()


def visualizer(mpc_code, bands, clear_dir):
    """
    Generates a series of histograms of monthly epochs across the various 
    years spanned by the MPC and WISE datasets of a given object.
    Arguments: mpc_file (str) -- the name of the MPC txt file
               wise_file (str) -- the name of the WISE tbl file
               clear_dir (bool) -- whether or not to clear the output folder
    Returns: a series of histograms in the output folder. Will also run
            the full pring out of the comparer function in 
            mpc_wise_functions.py
    """

    mpc_file, wise_file = return_input_files(mpc_code, bands)
    unique_epochs = comparer(mpc_code, bands, True)

    mpc_data = list(MPC_parser(mpc_code).keys())
    wise_data = list(WISE_parser(mpc_code, bands).keys())
    unique_data = list(unique_epochs.keys())

    if clear_dir:
        files = os.listdir("output_plots")
        for f in files:
            os.remove("output_plots/" + f)

    #Trimmed month list for MPC Data
    trimmed_month_mpc = []
    trimmed_year_mpc = []

    for i in mpc_data:
        trimmed_month_mpc.append(i[5:7])
    for i in mpc_data:
        trimmed_year_mpc.append(i[:4])

    #Trimmed month list for WISE Data
    trimmed_month_wise = []
    trimmed_year_wise = []

    for i in wise_data:
        trimmed_month_wise.append(i[5:7])
    for i in wise_data:
        trimmed_year_wise.append(i[:4])

    #Counting observations in each month for MPC
    mpc_count = dict(pd.Series(trimmed_month_mpc).value_counts())
    months = {"01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr", "05": "May", 
            "06": "Jun", "07": "Jul", "08": "Aug", "09": "Sep", "10": "Oct", 
            "11": "Nov", "12": "Dec"}
    print("MPC epoch counts by month:")
    for month in months:
        if month not in mpc_count.keys():
            mpc_count[month] = 0
        print(months[month] + ":", mpc_count[month])
    print("Total MPC epochs:", sum(mpc_count.values()))

    #Counting observations in each month for WISE
    wise_count = dict(pd.Series(trimmed_month_wise).value_counts())
    print("WISE epoch counts by month:")
    for month in months:
        if month not in wise_count.keys():
            wise_count[month] = 0
        print(months[month] + ":", wise_count[month])
    print("Total WISE epochs:", sum(wise_count.values()))

    # MPC month plots by year per month
    mpc_year_data = {}
    for datum in mpc_data:
        year, month, day = datum[:4], datum[5:7], datum[8:10]
        if year not in mpc_year_data.keys():
            mpc_year_data[year] = [datum]
        else:
            mpc_year_data[year].append(datum)

    for year in mpc_year_data:
        month_plot(mpc_year_data[year], "MPC", year)

    # WISE month plots by year per month
    wise_year_data = {}
    for datum in wise_data:
        year, month, day = datum[:4], datum[5:7], datum[8:10]
        if year not in wise_year_data.keys():
            wise_year_data[year] = [datum]
        else:
            wise_year_data[year].append(datum)

    for year in wise_year_data:
        month_plot(wise_year_data[year], "WISE", year)

    # New epoch month plots by year per month
    unique_year_data = {}
    for datum in unique_data:
        year, month, day = datum[:4], datum[5:7], datum[8:10]
        if year not in unique_year_data.keys():
            unique_year_data[year] = [datum]
        else:
            unique_year_data[year].append(datum)

    for year in unique_year_data:
        month_plot(unique_year_data[year], "Unique", year)

    month_plot(mpc_data, "MPC", "all years")
    month_plot(wise_data, "WISE", "all years")
    month_plot(unique_data, "Unique", "all years")