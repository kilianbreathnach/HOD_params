import sys
from glob import glob
import numpy as np
from dat.h5_funcs import h5_arr

survey_cap = "CMASS/NGC"
cosmo = "WMAP"
spheresfile = "dat/out/{0}/{1}/mocks/mock_srch_pts.hdf5".format(survey_cap, cosmo)

if len(sys.argv) != 2:
    print "usage: rdz2ascii.py <mock no.: int>"

digits = len(sys.argv[1])

i = 0
name = ""
while i < (4 - digits):
    name += "0"
    i += 1

name += sys.argv[1]

mock_rdzs = glob("dat/out/{0}/{1}/mocks/rdz_down/{2}.hdf5".format(survey_cap, cosmo, name))

for i in range(len(mock_rdzs)):

    mock_rdz = h5_arr(mock_rdzs[i], "radecz")

    np.savetxt("dat/out/{0}/{1}/mocks/rdz_down/ascii/{2}.dat".format(survey_cap, cosmo, name), mock_rdz)
