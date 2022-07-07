import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

data_object = Table.read("ne_inputs/161989_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

data_object2 = Table.read("mcmc_inputs/161989_2bands.tbl", format='ipac')
mjd2 = data_object2['mjd']
w12 = data_object2['w1flux']
w22 = data_object2['w2flux']

#plt.scatter(mjd2, w12, label='all w1')
plt.scatter(mjd2, w22, label='all w2')
#plt.scatter(mjd, w1, label='new w1')
plt.scatter(mjd, w2, label='new w2')
plt.legend()
plt.show()