from MPC_WISEcat_comparer import *
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd

mpc_data = '161989.txt'
wise_data = 'table_irsa_catalog_search_results-2.csv'
mpc_data, wise_data, new_epochs = comparison(mpc_data, wise_data)
print(new_epochs)


print(mpc_data)
print(wise_data)

#Trimmed month list for MPC Data
trimmed_month_mpc = []

for i in mpc_data:
    trimmed_month_mpc.append(i[5:7])

print(trimmed_month_mpc)

#Trimmed month list for WISE Data
trimmed_month_wise = []

for i in wise_data:
    trimmed_month_wise.append(i[5:7])

print(trimmed_month_wise)

#Counting observations in each month for MPC
mpc_count = pd.Series(trimmed_month_mpc).value_counts()
print("MPC Jan:" ,mpc_count["01"])
print("MPC Feb:" ,mpc_count["02"])
print("MPC March:" ,mpc_count["03"])
print("MPC April:" ,mpc_count["04"])
print("MPC May:" ,mpc_count["05"])
print("MPC June:" ,mpc_count["06"])
print("MPC July:" ,mpc_count["07"])
print("MPC Aug:" ,mpc_count["08"])
print("MPC Sept:" ,mpc_count["09"])
print("MPC Oct:" ,mpc_count["10"])
print("MPC Nov:" ,mpc_count["11"])
print("MPC Dec:" ,mpc_count["12"])

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



plt.hist(month_obs_mpc, alpha=0.5, label= 'MPC Data')
#plt.hist(y, alpha =0.5, label = 'WISE Data')
plt.xticks(np.arange(12), months)
plt.legend(loc='upper right')
plt.show()

#Count observations in each month

#for i in trimmed_month_mpc:
#     list = trimmed_month_mpc.count(i)
#  print(list)
