# WISEPY

WISEPY is a code base that serves to manipulate data from the WISE GATOR catalog and the International Astronomical Union's Minor Planet Center database. As a facet of the NEOWISE mission, it principally acts as a tool through which to recover missing observations from the WISE catalog through manual methods.

## General Requirements

WISEPY makes use of Linux commands when interacting with files, and as such, is not written to be run on a non-Linux based environment. Python is the only language used and I recommend using a conda installation to set up all of the packages used in this tool (Anaconda can be installed here: https://www.anaconda.com/). Apple computers likely already come with a version of Python installed. **DO NOT USE THIS!** Make sure your Python is up to date.

WISEPY also makes use of SAOImageDS9 through Linux commands for interacting with FITS files. Make sure DS9 is installed (it can be downloaded here: https://sites.google.com/cfa.harvard.edu/saoimageds9). In order for WISEPY to work, you must use the Linux executable build so that command line arguments can be used. For Mac users, do not use the application build of DS9.

In addition to the above points, make sure that your pathing is correct so DS9 can be used via the terminal. Additionally, if you are using VSCode, make sure that you are using the correct Python interpreter.

