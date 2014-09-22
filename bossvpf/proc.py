import os
import subprocess
import requests
from glob import glob
import numpy as np
from dat.cartesian_cosmo import mk_coords
from dat.h5_funcs import mk_h5, fits2h5
from dat.nbar_zrange import process_nbar


sdss_dir = "./dat/in/sdss3"
cmass_dir = "./dat/in/CMASS_DATA"

# CMASS
# - NGC

sdss_dir = "./dat/in/sdss3"
CMASS_dir = "./dat/in/CMASS_DATA"

out_dir = "./dat/out/CMASS/NGC"
cosmo = "WMAP"
cosmo_dir = "{0}/{1}".format(out_dir, cosmo)

# get the zrange around the mean nbar
process_nbar("{0}/nbar-cmass-dr11v1-N-Anderson.dat".format(sdss_dir),
             "{0}/{1}/nbar_zrange.json".format(out_dir, cosmo), cosmo,
             radeczfile="{0}/full_radecz.hdf5".format(out_dir))


