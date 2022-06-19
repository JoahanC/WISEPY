from mpc_wise_functions import *
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import datetime


def month_plot(dates, title):
    trimmed_dates = []
    for date in dates:
        year = date[:4]
        month = date[5:7]
        day = date[8:10]
        date_string = year + month + day
        trimmed_dates.append(date_string)

    dtm = lambda x: int(x[4:6])
    months = list(map(dtm, trimmed_dates))

    fig, ax = plt.subplots()
    bins = np.arange(1,14)
    ax.hist(months, bins = bins, edgecolor="k", align='left')
    ax.set_xticks(bins[:-1])
    ax.set_title(title)
    ax.set_xticklabels([datetime.date(1900,i,1).strftime('%b') for i in bins[:-1]] )
    plt.savefig(f"outputs/{title}")


mpc_data = '161989.txt'
wise_data = 'table_irsa_catalog_search_results.tbl'
unique_epochs = data_comparer(mpc_data, wise_data)
print(unique_epochs)

wise_data = list(WISE_parser(wise_data).keys())
mpc_data = list(MPC_parser(mpc_data).keys())
unique_data = list(unique_epochs.values())

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
months = ["01","02","03","04","05","06","07","08","09","10","11","12"]
print("MPC Jan:", mpc_count["01"])
print("MPC Feb:", mpc_count["02"])
print("MPC March:", mpc_count["03"])
print("MPC April:" ,mpc_count["04"])
print("MPC May:" ,mpc_count["05"])
print("MPC June:" ,mpc_count["06"])
print("MPC July:" ,mpc_count["07"])
print("MPC Aug:" ,mpc_count["08"])
print("MPC Sept:" ,mpc_count["09"])
print("MPC Oct:" ,mpc_count["10"])
print("MPC Nov:" ,mpc_count["11"])
print("MPC Dec:" ,mpc_count["12"])
print(mpc_count)
tot_count = []
#months = list(map(months, dates))

month_plot(mpc_data, "mpc")
month_plot(wise_data, "wise")
month_plot(unique_data, "unique")


#Counting observations in each month for WISE
wise_count = pd.Series(trimmed_month_wise).value_counts()
print("WISE Jan:" ,wise_count["01"])
print("WISE Feb:" ,wise_count["02"])

#number of observations per month into list for mpc
month_obs_mpc = [96, 481, 588, 146, 40, 26, 9, 37, 30, 15, 13, 25]
month_obs_wise = [10, 124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

x = month_obs_mpc
y = month_obs_wise


#labels = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#plt.bar(month_obs_mpc)
#plt.xticks(month_obs_mpc, labels)



"""
plt.hist(mpc_count, alpha=0.5, label= 'MPC Data')
#plt.hist(y, alpha =0.5, label = 'WISE Data')
plt.xticks(np.arange(12), months)
plt.legend(loc='upper right')
plt.show()

#Count observations in each month

#for i in trimmed_month_mpc:
#     list = trimmed_month_mpc.count(i)
#  print(list)"""
