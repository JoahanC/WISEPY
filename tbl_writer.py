from astropy.table import Table
from mpc_wise_functions import comparer


mpc_file = "161989.txt"
wise_file = "table_irsa_catalog_search_results.tbl"
new_epochs = comparer(mpc_file, wise_file, False)

source_id_list = []
for epoch in new_epochs:
    source_id_list.append(new_epochs[epoch][0])
name_list = ["cacus"] * len(source_id_list)
dat = Table([name_list, source_id_list], names=["name", "source_id"])
dat.write("table.tbl", format="ipac")