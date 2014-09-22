import sys
from glob import glob
import numpy as np
from nbar_utils import nzobj


def process_mocks(dattop, sample, cap, cosmology, cutoffval):

    # First get the polynomial fit to the average of all the mock number densities
#    mock_avg_poly()

    # initialise an nzobj
    mock_nz = nzobj(dattop, sample, cap, cosmology, cutoffval)

    # load the json data into the nzobj
    mock_nz.read_json()

    # give the nzobj the polynomial coefficients to use for downsampling
    mock_nz.popt = np.loadtxt(dattop + "/out/{0}/{1}/{2}/mocks/avg_poly_coeffs.dat".format(sample, cap, cosmology))

    # glob all the mocks
    mockfns = glob("{0}/in/boss/{1}/{2}/mocks/*.rdz".
                                 format(dattop, sample, cap))

    # now process all those bad boys
    for mockfn in mockfns:

        mockN = mockfn.split("/")[-1].split("_")[1].split(".")[0]

        print "\n", "mock\t", mockN, "\n"
        mock_nz.radecz = np.loadtxt(mockfn)

        mock_nz.rdz_down("/mocks/rdz/{0}.dat".format(mockN))
