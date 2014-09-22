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


def plot_nz(fn):

    werds = fn.split("/")
    outi = werds.index("out")

    dattop = "/".join(werds[:outi])
    sample = werds[outi + 1]
    cap = werds[outi + 2]
    cosmo = werds[outi + 3]

    nob = nzobj(dattop, sample, cap, cosmo, 0.65)

    nob.read_json()

    nob.radecz = np.loadtxt(fn)

    nob._compute_z_nbar()

    fig = plt.figure()

    ax = fig.add_subplot(111)

    l, = ax.plot(nob.z_cen, nob.nbar)

    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$ n(z)$")

    picname = ''.join(os.path.basename(fn).split(".")[:-1]) + ".pdf"
    plotaddress = os.path.dirname(fn) + "/plots/" + picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument"

    plot_nz(sys.argv[1])
