import os
import sys
from glob import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

from astropy.io import fits
from nbar_utils import nzobj


def plot_nz(datrnk, fndown):

    sample = "CMASS"
    cap = "NGC"

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot corrected
    nob = nzobj(datrnk, sample, cap, "WMAP", 0.65)

    nob.get_zrange()

    ax.scatter(nob.nbar_corr[:, 0], nob.nbar_full,
               color='m', label=r"$n_{corrected}$")

    # plot raw data
    galfits = fits.open(glob("{0}/in/boss/{1}/{2}/*.fits".
                                 format(datrnk, sample, cap))[0])

    # array of all the RAs Decs and redshifts
    radecz = np.dstack((galfits[1].data["RA"],
                        galfits[1].data["DEC"],
                        galfits[1].data["Z"]))[0]

    buf = len(nob.zcen)

    # give the rdz to the nzobj
    nob.radecz = radecz

    nob._compute_z_nbar(buf=buf)

    ldat, = ax.plot(nob.z_cen, nob.nbar, label=r"$n_{data}$")

    # plot downsampled nbar
    nob.radecz = np.loadtxt(fndown)

    nob._compute_z_nbar(buf=5)

    nob.nbar[:5] = 0.
    nob.nbar[-5:] = 0.

    ldown, = ax.plot(nob.z_cen, nob.nbar,
                     color='r', label=r"$n_{downsampled}$")

    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$ n(z)$ $(Mpc / h)^{-3}$")
    ax.set_title("CMASS NGC Number Densities")

    ax.set_xlim([0.3, 0.75])
    ax.set_ylim([0., 0.0005])

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    picname = "nbar_talk.pdf"
    plotaddress = picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))



if __name__ == "__main__":

    plot_nz(sys.argv[1], sys.argv[2])
