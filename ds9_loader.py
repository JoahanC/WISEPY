import os 



ds9_script = "ds9 -tile "
files = ["56769b152-w1-int-1b_ra218.094537_dec-16.22159_asec600.000.fits",
         "56769b152-w2-int-1b_ra218.094537_dec-16.22159_asec600.000.fits",
         "56816a125-w1-int-1b_ra220.016094_dec-15.68386_asec600.000.fits",
         "56816a125-w2-int-1b_ra220.016094_dec-15.68386_asec600.000.fits",
         "56844a124-w1-int-1b_ra221.20409_dec-15.33466_asec600.000.fits",
         "56844a124-w2-int-1b_ra221.20409_dec-15.33466_asec600.000.fits",
         "56848a123-w1-int-1b_ra221.374803_dec-15.28348_asec600.000.fits",
         "56848a123-w2-int-1b_ra221.374803_dec-15.28348_asec600.000.fits",
         "56865b147-w1-int-1b_ra222.146896_dec-15.04868_asec600.000.fits",
         "56865b147-w2-int-1b_ra222.146896_dec-15.04868_asec600.000.fits",
         "56884a122-w1-int-1b_ra222.924815_dec-14.80676_asec600.000.fits",
         "56884a122-w2-int-1b_ra222.924815_dec-14.80676_asec600.000.fits",
         "56964a119-w1-int-1b_ra226.44876_dec-13.64611_asec600.000.fits",
         "56964a119-w2-int-1b_ra226.44876_dec-13.64611_asec600.000.fits",
         "57000a117-w1-int-1b_ra228.069343_dec-13.07775_asec600.000.fits",
         "57000a117-w2-int-1b_ra228.069343_dec-13.07775_asec600.000.fits"]

for file in files:
    sid, w_band = file[:9], file[10:12]
    ds9_script += "wise_images/161989/" + file + ' -regions '
    ds9_script += f"regions/{sid}_{w_band}.reg "
os.popen(ds9_script)
