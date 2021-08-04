#Data Handling
import numpy as np
import itertools
import argparse 
import PhysicalConstantsCGS as const    
import pickle
from scipy.optimize import curve_fit

#System modules
import time
import os
import sys
import warnings
import shutil
import subprocess

#Astropy Data Handling
from astropy import units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
import scipy

#Visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmasher as cmr
mpl.style.use('classic')
# To reproduce plotting settings Menon used
mpl.rc_file('./matplotlibrc',use_default_template=False)


##### Some global quantities #########
arcsec_to_degree = 1./3600.
list_of_galaxies = ['NGC_0628','NGC_1313','NGC_1566','NGC_3344','NGC_3627',
        'NGC_3738','NGC_4449','NGC_5194',
        'NGC_5253','NGC_5457','NGC_6503','NGC_7793']
       

mosaic_links = \
['https://archive.stsci.edu/hlsps/legus/mosaics/ngc628/hlsp_le\
gus_hst_acs_ngc628-mosaic_f435w_v1_sci.fits',
'https://archive.stsci.edu/hlsps/legus/mosaics/ngc1313/hlsp_\
legus_hst_acs_ngc1313-mosaic_f435w_v1_sci.fits', #435,555,814
'https://archive.stsci.edu/hlsps/legus/ngc1566/ngc1566_drc.tar.gz',#any
'https://archive.stsci.edu/hlsps/legus/ngc3344/ngc3344_drc.tar.gz',#any
'https://archive.stsci.edu/hlsps/legus/ngc3627/ngc3627_drc.tar.gz',#any
'https://archive.stsci.edu/hlsps/legus/ngc3738/ngc3738_drc.tar.gz',#606,814
'https://archive.stsci.edu/hlsps/legus/ngc4449/ngc4449_drc.tar.gz',#435,555,814
'https://archive.stsci.edu/hlsps/legus/mosaics/ngc5194_ngc5195/h\
lsp_legus_hst_acs_ngc5194-ngc5195-mosaic_f435w_v1_sci.fits',
'https://archive.stsci.edu/hlsps/legus/ngc5253/ngc5253_drc.tar.gz',#435,555,814
'https://archive.stsci.edu/hlsps/legus/mosaics/ngc5457/hlsp_legus\
_hst_wfc3_ngc5457-mosaic_f336w_v1_sci.fits',
'https://archive.stsci.edu/hlsps/legus/ngc6503/ngc6503_drc.tar.gz',#any
'https://archive.stsci.edu/hlsps/legus/mosaics/ngc7793/hlsp_legus_\
hst_acs-wfc3_ngc7793-mosaic_f555w_v1_sci.fits'
]

# Pickle Data Handling
def saveObj(obj, name):
    """
    Save a pickle object.

    INPUTS:
    ----------
    obj      - the name of the data object to be saved
    name     - the name of the pickle object

    """

    os.system("touch " + name + ".pkl")
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(obj, f, protocol=2)


def loadObj(name):
    """
    Load a pickle object.

    INPUTS:
    ----------
    name     - the name of the pickle object

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
