import os
import sys
import subprocess
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


def plot_vpfs(vpf1, vpfdat, daterr, picname):

    # Get file data
    vpf1 = np.loadtxt(vpf1)
#    vpf2 = np.loadtxt(vpf2)
    vpfdat = np.loadtxt(vpfdat)
    daterr = np.loadtxt(daterr)

    fig = plt.figure()
    ax = fig.add_subplot(111)


    l1, = ax.plot(vpf1[:, 0],
                  np.log10(np.abs(vpf1[:, 1] - vpfdat[:, 1]) / daterr[:, 1]),
                  color='b')
#    l2, = ax.plot(vpf2[:, 0],
#                  np.log10(np.abs(vpf2[:, 1] - vpfdat[:, 1]) / daterr[:, 1]),
#                  color='g')


    # annotate

    ax.set_ylabel(r"$log(\frac{|\Delta P_{0}(r)|}{\sigma_{P_{0}}})$")
    ax.set_xlabel(r"$r$ $(Mpc)$")
#    ax.set_ylim([-2., 2.])

    plotaddress = picname

    plt.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))


if __name__ == "__main__":

    if len(sys.argv) != 6:
        print "usage: python vpf_sigma_diff.py <vpf1> <vpf2> \
              <vpf data> <data errors> <png output filename>"

    plot_vpfs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

