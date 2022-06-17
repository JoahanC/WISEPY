from MPC_WISEcat_comparer import *


mpc_data = '161989.txt'
wise_data = 'table_irsa_catalog_search_results-2.csv'
mpc_data, wise_data, new_epochs = comparison(mpc_data, wise_data)
print(new_epochs)