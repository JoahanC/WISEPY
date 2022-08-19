"""
This file is a basic tutorial to using WISEPY.
"""
from asteroid import Asteroid


"""
WISEPY is a tool that can be used to look at asteroid image data from the WISE
catalog and perform data analysis as well as cross referencing with the 
corresponding database for the given asteroid from the Minor Planet Center.
In order for data to be properly parsed by WISEPY, all input data must be stored
within the database_files folder. In the database_files folder, there are two
subfolders corresponding to the mpc and wise data. A specific naming scheme
must be followed. For MPC data, WISEPY will only recognized .txt files, and they
must be named in the following convention: [[packed_name]].txt. For WISE data,
WISEPY also follows a strict nameing scheme, only accepting files formatted
in the following manner: [[packed_name]]_[[band]].tbl. In this case, the [[band]]
corresponds to the amount of wavelength bands present in the .tbl file, which will
vary depending on which year of NEOWISE data you are looking at.
"""

"""
First we will initialize the Asteroid class which hosts all of WISEPYs functionality.
To initialize an Asteroid object, two arguments are needed. Firstly, the packed
name consistent with all of the input files in /database_files/, and secondly, a list
of all the different NEOWISE phases from which data was pulled. This list can be
constructed by looking at the wise input files for the chosen object and making a 
list of the bands correspoding to each object. For example: 02212_2band.tbl, 
02100_4band.tbl would be initialized as:
"""
cacus = Asteroid("02212", [2, 4])

"""
First we'll run the comparison feature WISEPY possesses. To do this, simply call the
run_comparison method. This will internally populate the asteroid object with all of
the unique observations found in the WISE dataset relative to the MPC dataset.
Additionally, WISEPY will generate SNR and flux plots for all detections and solely
unique detections in the /plots/ folder.
"""
cacus.run_comparison()
cacus.display_unique_observations()

"""
Unique detections can be viewed through WISEPYs DS9 interface. After making sure DS9
is installed, you can view either all of the unique detections images, or only a subset
defined in /loader_data/ slash using the ds9_viewer() or ds9_loader([[band]]),
respectively where the [[band]] corresponds to the NEOWISE image subset. ds9_loader
will only work if an existing .txt file corresponding to the band desired exists, even
if it is only blank.
"""
cacus.ds9_loader(2)
cacus.ds9_viewer()

"""
Note that the ds9_viewer() methods will allow you to view as many FITS images as you
want. However, you will be capped by the limitations of your computer and may crash
DS9 should you load to many. For this reason, you will be prompted to specify a range
of images to view, should you want to view only batches of DS9 code.
"""

"""
/loader_data/ hosts all of the detections can may be deemed valid after manual inspection
of FITS images. These detections are simply written to a .txt file manually and only
contain source_id stubs, such as 021492a121, 31513b124, etc. 
"""

"""
Individual dictionaries of these new detections can be retrieved with the get_* methods
of the Asteroid class. 
"""
cacus.get_mpc_observations(2)
cacus.get_wise_observations(2)

"""
Additionally, astropy tables can be written with either
only unique detections, or all valid detections. Tables written for new_detections will
be blank unless a specific set of source_id stubs are written to a .txt file in 
/loader_data/. NOTE that .txt files for all appropriate bands must be writen (even if 
there are no unique detections for that band) in /loader_data/ before using these methods, 
otherwise WISEPY will look for a file that doesn't exist. NOTE generate_new_table must be 
run before generate_full_table. Tables are outputted to the /observations/ folder in 
either `all` or `new`.
"""
cacus.generate_new_table()
cacus.generate_full_table()



