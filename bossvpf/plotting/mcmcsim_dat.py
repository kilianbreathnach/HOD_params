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

    fig = plt.figure()
    ax = fig.add_subplot(111)

    simfn = dattop + "/out/CMASS/NGC/WMAP/vpf/mcmc/new/frac_var.dat"

    vpf = np.loadtxt(simfn)
    zs = np.append(0., vpf[:, 0])
    meanvpf = np.append(1., vpf[:, 1])
    sigma = np.append(0., np.sqrt(vpf[:, 2]))

    l1, = ax.plot(zs, meanvpf,
                  color='g', ls='--', alpha=0.6,
                  label="simulation")

    ax.fill_between(zs, meanvpf + sigma, meanvpf - sigma,
                    facecolor='green', alpha=0.4)

    datfn = dattop + "/out/CMASS/NGC/WMAP/vpf/comp_0.05.vpf"
    vpf = np.loadtxt(datfn)

#    l2, = ax.plot(np.append(0., vpf[:, 0]),
#                  np.append(1., vpf[:, 1]),
#                  marker='o', color='b',
#                  alpha=0.5, label='data')

    errors = np.loadtxt(dattop + "/out/CMASS/NGC/WMAP/vpf/mocks/frac_var.dat", usecols=(1,))

    ax.errorbar(np.append(0., vpf[:, 0]),
                np.append(1., vpf[:, 1]),
                yerr=np.append(0., vpf[:, 1] * np.sqrt(errors)),
                fmt='.', lw=1, capsize=5, mew=1,
                label="data")
#                linestyle='none', marker='.', fmt='')
#                 marker='o', color='k', ecolor='k',
#                 markerfacecolor='b', label='data')
     #            linestyle=None) # fmt=None, color='b')

    ax.set_yscale('log', nonposy='clip')
    ax.set_xlim([0., 45.])

    ax.set_xlabel(r"$r$ $(h^{-1} Mpc)$")
    ax.set_ylabel(r"$P_{0}(r)$")
    ax.set_title("VPF for CMASS NGC Galaxies")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    picname = "mcmcsimandat.pdf"
    plotaddress = os.path.dirname(datfn) + "/plots/" + picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument: the address of the data directory"


    plot_vpfs(sys.argv[1])
