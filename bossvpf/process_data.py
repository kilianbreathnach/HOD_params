import os
from glob import glob
import numpy as np
from astropy.io import fits
from nbar_utils import nzobj


def process_data(dattop, sample, cap, cosmology, cutoffval):

    data_nz = nzobj(dattop, sample, cap, cosmology, cutoffval)

    data_nz.get_zrange()

    # get the galaxy data to compute number density from
    galfits = fits.open(glob("{0}/in/boss/{1}/{2}/*.fits".
                                 format(dattop, sample, cap))[0])

    # array of all the RAs Decs and redshifts
    radecz = np.dstack((galfits[1].data["RA"],
                        galfits[1].data["DEC"],
                        galfits[1].data["Z"]))[0]

    # give the rdz to the nzobj
    data_nz.radecz = radecz

    data_nz.fit_polynom()

    data_nz.rdz_down("/radecz.dat")

    data_nz.write_json()
