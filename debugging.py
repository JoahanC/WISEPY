from MPC_WISEcat_comparer import *


mpc_file = '161989.txt'
wise_file = 'table_irsa_catalog_search_results.tbl'

unique_epochs = data_comparer(mpc_file, wise_file)
