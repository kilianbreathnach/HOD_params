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


def plot_vpfs(dattop):

    mockfns = glob(dattop + "/out/CMASS/NGC/WMAP/vpf/mocks/*.vpf")

    fig = plt.figure()
    ax = fig.add_subplot(111)

    lines = []

    for fn in mockfns:

        vpf = np.loadtxt(fn)
        try:
            lines.append(ax.plot(np.append(0., vpf[:, 0]),
                                 np.append(1., vpf[:, 1]),
                                 color='b', alpha=0.05))
        except:
            continue

    simfn = dattop + "/out/CMASS/NGC/WMAP/vpf/sim_com_0.95.vpf"

    vpf = np.loadtxt(simfn)
    l1, = ax.plot(np.append(0., vpf[:, 0]),
                  np.append(1., vpf[:, 1]),
                  color='g', alpha=0.5)

    datfn = dattop + "/out/CMASS/NGC/WMAP/vpf/comp_0.05.vpf"

    errors = np.loadtxt(dattop + "/out/CMASS/NGC/WMAP/vpf/mocks/frac_var.dat", usecols=(1,))

    ax.errorbar(np.append(0., vpf[:, 0]),
                np.append(1., vpf[:, 1]),
                yerr=np.append(0., errors), fmt='o')

    ax.set_yscale('log')

    ax.set_xlabel(r"$r$ $(Mpc)$")
    ax.set_ylabel(r"$P_{0}(r)$")

    picname = "simandat.pdf"
    plotaddress = os.path.dirname(datfn) + "/plots/" + picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument: the address of the data directory"


    plot_vpfs(sys.argv[1])
