import os
import json
import numpy as np
from glob import glob
from astropy.io import fits
from astropy.coordinates import Distance
from astropy.cosmology import Planck13, WMAP5
from h5_funcs import arr2h5, h5_arr


"""
This file contains functions to take the BOSS data file for corrected number
density as a function of redshift and process it for our purposes.

process_nbar decides on a width around the redshift of maximum number density
of each survey that we focus on for our analysis (decided by Q below) and
computes the overall average in that redshift range. It stores
information to a json output file and then uses this information to process
either the data or a mock data file.
"""


def calc_shell_vol(dis_func, z_near, z_far, zcen):
    """
    Computes the volume of a spherical shell of at redshift zcen, with edges
    at redshifts z_near and z_far, given a cosmological distance conversion
    function. Results are in (h^-1 Mpc)^3.
    """

    r = dis_func(zcen).value
    r_near = dis_func(z_near).value
    r_far = dis_func(z_far).value
    dr = r_far - r_near

    return 4 * np.pi * (r ** 2) * dr


def process_nbar(nbarfile, nz_dict_file, cosmology, radeczfile=None):
    """
    Parameters
    ---------

    nbarfile : str
        the path to and name of the corrected nbar file
    nz_dict_file : str
        path to and name of the json file with the nbar dict
    cosmology : str, "WMAP" or "Planck"
        the cosmology to compute shell volumes with
    radeczfile : str, "data" or "mock"
        the data or mock file to process
    """

    # magic number for width around maximum
    Q = 0.65
    # magic number for shell vol computation
    Nfrac = (6769.0358 * np.pi) / 129600

    if cosmology == "Planck":
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif cosmology == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5
    comv = cosmo.comoving_distance

    nbar_corr = np.loadtxt(nbarfile)
    nz_dict = {"tophat height for zrange": Q}

    # Cut out the first bit of crap (works for CMASS, dunno about LOWZ)
    ind03 = np.abs(nbar_corr[:, 0] - 0.3).argmin()

    nbar_corr = nbar_corr[ind03:, :]

    zcen = nbar_corr[:, 0]
    z_near = nbar_corr[:, 1]
    z_far = nbar_corr[:, 2]
    corr_gal_counts = nbar_corr[:, 6]

    nbar = []
    shell_vols = []

    for i in range(len(zcen)):

        shell_vols.append(Nfrac * calc_shell_vol(comv, z_near[i], z_far[i], zcen[i]))
        nbar.append(corr_gal_counts[i] / shell_vols[i])

    nbar = np.array(nbar)
    shell_vols = np.array(shell_vols)

    # Find nbar peak and index
    max_nbar = np.max(nbar)
    max_i = int(np.where(nbar == max_nbar)[0])

    nz_dict["max_nbar_corr"] = max_nbar
    nz_dict["nbar_corr_tophat"] = Q * max_nbar
    nz_dict["z_nbar_max"] = zcen[max_i]

    # get the interval edge indices
    L = np.abs(nbar[:max_i] - max_nbar * Q).argmin()
    R = max_i + np.abs(nbar[max_i:] - max_nbar * Q).argmin()

    nbar = nbar[L:R + 1]
    shell_vols = shell_vols[L:R + 1]

    nz_dict["zlo"] = zcen[L]
    nz_dict["zhi"] = zcen[R]

    nz_dict["avg_nbar_corr"] = np.average(nbar)
    nz_dict["total_shell_vol"] = np.sum(shell_vols)

    if radeczfile:

        radecz = h5_arr(radeczfile, "radecz")

        # Make the redshift cut in the nbar array with right cosmology
        nbar_corr = nbar_corr[(nz_dict["zlo"] <= nbar_corr[:, 0]) * \
                            (nbar_corr[:, 0] <= nz_dict["zhi"])]

        # Get binning those observed galaxies
        zbinedges = np.append(nbar_corr[0, 1], nbar_corr[:, 2])

        # Find the counts per bin and convert to nbar
        H = np.histogram(radecz[:, 2], bins=zbinedges)
        hist_nbar = H[0] / shell_vols

        if not radeczfile.split('/')[-2] == "mocks_hierarchical":
            # save the average downsampled value if it's the data file
            nz_dict["avg_nbar_down"] = np.average(hist_nbar)

            # The number to downsample to in each bin
            # (multiply bin number by the relative fraction determined from
            #  corrected distribution of nbar)
            nz_dict["nbar_data_tophat"] = 0.95 * nz_dict["nbar_corr_tophat"] * (nz_dict["avg_nbar_down"] / nz_dict["avg_nbar_corr"])
            factor_arr = nz_dict["nbar_data_tophat"] / hist_nbar

        # if we are dealing with a mockfile, then there is an extra factor to make
        # the average equal that of the data
        elif radeczfile.split('/')[-2] == "mocks_hierarchical":

            # have to open the existing json file
            jf = open(nz_dict_file)
            nz_dict = json.load(jf)

            factor_arr = nz_dict["nbar_data_tophat"] / hist_nbar

            jf.close()

        num_down = np.rint(factor_arr * H[0])
        num_down = num_down.astype(int)

        # make a mask for the final array for analysis within the redshift limits
        finmask = np.array(radecz.shape[0] * [False])

        for i, nd in enumerate(num_down):
            """Turn on the right amount of galaxies in each bin."""
            zbin_ids = np.where(((zbinedges[i] < radecz[:, 2]) * (radecz[:, 2] <= zbinedges[i + 1])) == True)

            if zbin_ids[0].shape[0] == 0:
                continue

            keep = np.random.choice(zbin_ids[0], size=nd, replace=False)

            finmask[keep] = True

        radecz = radecz[finmask]

        if not radeczfile.split('/')[-2] == "mocks_hierarchical":
            # and save downsampled array to a hdf5 file
            arr2h5(radecz, "{0}/radecz_down.hdf5".format(os.path.dirname(nz_dict_file)), "radecz", mode='w')

        elif radeczfile.split('/')[-2] == "mocks_hierarchical":
            # save mocks to a hdf5 file
            mock_no = radeczfile.split('/')[-1].split('.')[0]
            arr2h5(radecz, "{0}/mocks/rdz_down/{1}.hdf5".format(os.path.dirname(nz_dict_file), mock_no), "radecz", mode='w')


    if not radeczfile.split('/')[-2] == "mocks_hierarchical":
        # don't save the json if we're working on a mock
        nf = open(nz_dict_file, 'w')
    
        json.dump(nz_dict, nf, sort_keys=True, indent=4, separators=(',', ':\t'))
    
        nf.close()
