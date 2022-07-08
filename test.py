import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

data_object = Table.read("new_inputs/2002_QF15_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

plt.scatter(mjd, w2, label='W2 band')
plt.ylabel('Flux')
plt.xlabel('MJD')
plt.legend()
plt.show()