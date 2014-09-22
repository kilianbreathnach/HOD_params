import os
import sys
from glob import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

from cosmo_stuff import get_inv_efunc, get_comv


inv_efunc = get_inv_efunc("WMAP")
comv = get_comv("WMAP")

def plot_vpfs(dattop):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    z = 0.5342
    fibrad = np.radians(62. / 3600) * comv(z).value

    velfac = inv_efunc(z) * (1 + z)

    dat = np.loadtxt(dattop + "/in/boss/CMASS/NGC/cmass-dr11v2-N-Anderson-disp_los_pm.dat")

    mock = np.loadtxt(dattop + "/out/sim/simdown.dat")
    mockxy = mock[:, :2]
    gal_baum = cKDTree(mockxy)

    collisions = gal_baum.query_ball_tree(gal_baum, fibrad)

    dlos = []
    redshifty = []

    for i, collision in enumerate(collisions):

        if len(collision) > 1:

            for gal in collision:

                if not gal == i:

                    dlos.append(np.abs(mock[i, 2] - mock[gal, 2]))
                    rsd1 = mock[i, 2] + 0.01 * mock[i, 5] * velfac
                    rsd2 = mock[gal, 2] + 0.01 * mock[gal, 5] * velfac

                    redshifty.append(np.abs(rsd1 - rsd2))


    dlos = np.array(dlos)
    redshifty = np.array(redshifty)

    binedges = np.append(np.arange(0., 40, 1),
                         np.arange(40., redshifty.max() + 5, 5))

    dlosH = np.histogram(dlos, bins=binedges)
    redsH = np.histogram(redshifty, bins=binedges)
    datsH = np.histogram(dat, bins=binedges)

    r = binedges[:-1] + 2.5

    ldlos, = ax.plot(r, 0.5 * dlosH[0], label=r"mock $d_{LOS}$")
    lreds, = ax.plot(r, 0.5 * redsH[0], label=r"mock $d_{LOS}$ with RSD")
    ldats, = ax.plot(r, 0.5 * datsH[0], label=r"data $d_{LOS}$")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    ax.set_xlabel(r"$r$ $(h^{-1} Mpc)$")

    picname = "fibre_dlos_full.pdf"
    plotaddress = picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument: the address of the data directory"


    plot_vpfs(sys.argv[1])
