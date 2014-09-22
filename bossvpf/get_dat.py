import os
import subprocess
import request

"""
This script pulls all the data from the boss


"""

def get_dat():

    print "\nDownloading data...\n\n"

    # Location of the data on the interwebz
    broiler_dat = "http://broiler.astrometry.net/~kilian/tinker/proj_dat/"
    sdss_site = "http://data.sdss3.org/sas/bosswork/boss/lss/"

    # Prompts for the sdss3 data username and password,
    # please contact the authors for information on these
    username = raw_input("please enter the username for the BOSS\
                                        data page: ")
    password = raw_input("and the corresponding password: ")


    # CMASS data
    # - NGC data
    ############

    # File names of all the input data
    nbar = "nbar-cmass-dr11v1-N-Anderson.dat"
    fitsfile = "cmass-dr11v1-N-Anderson.dat.fits"
    rands = "randoms_1E6_cmass_N_dr11v1.dat"

    # Create necessary directory structure if not in place
    sdss_dir = "./dat/in/sdss3"
    if not os.path.exists(sdss_dir):
        os.makedirs(sdss_dir)

    cmass_dir = "./dat/in/CMASS_DATA"
    if not os.path.exists(cmass_dir):
        os.makedirs(cmass_dir)

    # First, get the corrected number densities
    r_nbar = requests.get(sdss_site + nbar, auth=(username, password))

    # move them to the sdss3 data directory
    f = open("{0}/{1}".format(sdss_dir, nbar), 'w')
    f.write(r_nbar.text)
    f.close()

    # get the ra, dec and redshift information from sdss fits file
    r_radecz = requests.get(sdss_site + fitsfile, auth=(username, password))

    f = open("{0}/{1}".format(sdss_dir, fitsfile), 'w')
    f.write(r_radecz.text)
    f.close()

    # get the random search points and move them to CMASS_DATA
    r_rands = requests.get(broiler_dat + rands)

    f = open("{0}/{1}".format(cmass_dir, rands), 'w')
    f.write(r_rands.text)
    f.close()

    # get the bad points
    r_bad1 = requests.get(broiler_dat + "north_block_badfield.dat")
    r_bad2 = requests.get(broiler_dat + "north_block_brightstar.dat")
    r_bad3 = requests.get(broiler_dat + "north_block_outside.dat")

    f = open("{0}/{1}".format(cmass_dir, "north_block_badfield.dat"), 'w')
    f.write(r_bad1.text)
    f.close()

    f = open("{0}/{1}".format(cmass_dir, "north_block_brightstar.dat"), 'w')
    f.write(r_bad2.text)
    f.close()

    f = open("{0}/{1}".format(cmass_dir, "north_block_outside.dat"), 'w')
    f.write(r_bad3.text)
    f.close()


