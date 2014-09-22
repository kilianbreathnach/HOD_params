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

from nbar_utils import nzobj


def plot_nz(dirnam):

    werds = dirnam.split("/")
    outi = werds.index("out")

    dattop = "/".join(werds[:outi])
    sample = werds[outi + 1]
    cap = werds[outi + 2]
    cosmo = werds[outi + 3]

    allmocks = glob(dirnam + "/*.dat")

    nob = nzobj(dattop, sample, cap, cosmo, 0.65)
    nob.read_json()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    lines = []
    i = 0 # just get the first 100

    for fn in allmocks:

        i += 1
        if i > 100:
            break

        nob.radecz = np.loadtxt(fn)
        nob._compute_z_nbar()

        lines.append(ax.plot(nob.z_cen, nob.nbar, color='b', alpha=0.05))

    datfn = dattop + "/".join([sample, cap, cosmo]) + "/radecz.dat"

    nob.radecz = np.loadtxt(datfn)
    nob._compute_z_nbar()

    lines.append(ax.plot(nob.z_cen, nob.nbar, color='r', alpha=0.5))

    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$ n(z)$")
    ax.set_ylim([0.0002, 0.0003])

    picname = "mock_dat_nbars.pdf"
    plotaddress = dirnam + "/plots/" + picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument"

    plot_nz(sys.argv[1])
