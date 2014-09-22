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


def plot_allvpfs(dattop):

    mockfns = glob(dattop + "/out/CMASS/NGC/WMAP/vpf/mocks/*.vpf")

    fig = plt.figure()
    ax = fig.add_subplot(111)

    lines = []

    i = 0

    for fn in mockfns:

        vpf = np.loadtxt(fn)

        if i == 0:
            lines.append(ax.plot(vpf[:, 0], vpf[:, 1],
                         color='b', alpha=0.01, label="mocks"))
        i += 1

        try:
            lines.append(ax.plot(vpf[:, 0], vpf[:, 1],
                         color='b', alpha=0.05))
        except:
            continue

    datfn = dattop + "/out/CMASS/NGC/WMAP/vpf/comp_0.05.vpf"

    vpf = np.loadtxt(datfn)
    lines.append(ax.plot(vpf[:, 0], vpf[:, 1],
                 color='r', alpha=0.5, label="data"))


    ax.set_yscale('log')
    ax.set_xlim([0., 45.])

    ax.set_xlabel(r"$r$ $(h^{-1} Mpc)$")
    ax.set_ylabel(r"$P_{0}(r)$")
    ax.set_title("Mock VPF comparison")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    picname = "mockandat.pdf"
    plotaddress = os.path.dirname(datfn) + "/plots/" + picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument: the address of the data directory"


    plot_allvpfs(sys.argv[1])
