from MPC_WISEcat_comparer import *
import matplotlib.pyplot as plt

mpc_data = '161989.txt'
wise_data = 'table_irsa_catalog_search_results-2.csv'
mpc_data, wise_data, new_epochs = comparison(mpc_data, wise_data)
print(new_epochs)


print(mpc_data)
print(wise_data)


mpc_data_list = mpc_data.tolist()
wise_data_list = wise_data.tolist()

#Trimmed month list for MPC Data
trimmed_month_mpc = []

for i in mpc_data_list:
    trimmed_month_mpc.append(i[5:7])

print(trimmed_month_mpc)

#Trimmed month list for WISE Data
trimmed_month_wise = []

for i in wise_data_list:
    trimmed_month_wise.append(i[5:7])

print(trimmed_month_wise)