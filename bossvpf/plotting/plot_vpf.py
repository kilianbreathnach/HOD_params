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


def plot_vpf(fn):

    vpf = np.loadtxt(fn)

    fig = plt.figure()

    ax = fig.add_subplot(111)

    l, = ax.plot(vpf[:, 0], vpf[:, 1])

    ax.set_yscale('log')

    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$P_{0}(r)$")

    picname = ''.join(os.path.basename(fn).split(".")[:-1]) + ".pdf"
    plotaddress = os.path.dirname(fn) + "/plots/" + picname

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument"

    plot_vpf(sys.argv[1])
