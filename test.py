import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

data_object = Table.read("ne_inputs/23606_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

plt.plot(mjd, w1)
plt.plot(mjd, w2)
plt.savefig('as1.png')
plt.clf()

data_object = Table.read("ne_inputs/161989_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

plt.plot(mjd, w1)
plt.plot(mjd, w2)
plt.savefig('cacus.png')
plt.clf()

data_object = Table.read("ne_inputs/5693_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

plt.plot(mjd, w1)
plt.plot(mjd, w2)
plt.savefig('5693.png')
plt.clf()

data_object = Table.read("ne_inputs/2002_QF15_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

plt.plot(mjd, w1)
plt.plot(mjd, w2)
plt.savefig('2002.png')
plt.clf()

data_object = Table.read("ne_inputs/7335_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

plt.plot(mjd, w1)
plt.plot(mjd, w2)
plt.savefig('7335.png')
plt.clf()

data_object = Table.read("ne_inputs/85713_2bands.tbl", format='ipac')
mjd = data_object['mjd']
w1 = data_object['w1flux']
w2 = data_object['w2flux']

plt.plot(mjd, w1)
plt.plot(mjd, w2)
plt.savefig('85713.png')
plt.clf()
