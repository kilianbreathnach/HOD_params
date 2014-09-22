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


def plot_vpfs(wpfile, wpcovar, wpsim, wpsimdown, wpsimfibre, wpsimfibre2,picname):

    # Get date to add to plot
    process = subprocess.Popen(['date', '+%d-%b-%Y'], stdout=subprocess.PIPE)
    date, err = process.communicate()
    date = date.split()[0].split("-")

    # Get file data
    wpfile = np.loadtxt(wpfile)
    wpcovar =np.loadtxt(wpcovar)
    wpsim = np.loadtxt(wpsim)
    wpsimdown = np.loadtxt(wpsimdown)
    wpsimfibre = np.loadtxt(wpsimfibre)
    wpsimfibre2 = np.loadtxt(wpsimfibre2)

#    clustering = np.loadtxt(simclust, usecols=(0, 4))


    # Set up axis dimensions
    wp_ax = plt.axes([0.1, 0.35, 0.9, 0.65])
    res_ax = plt.axes([0.1, 0.1, 0.9, 0.2])

    rows = wpfile.shape[0]

    # plot wp from wp_covar from real mock catalogue against .clustering values

    errors = np.log10(np.diag(wpcovar.reshape((rows, rows))))

    fit1, = wp_ax.plot(wpsim[:, 0], wpsim[:, 1],
                       color='g', label=r"$M_{sat} = 0.5 \times M1$")
    fit2, = wp_ax.plot(wpsimdown[:, 0], wpsimdown[:, 1],
                       color='m', label=r"$M_{sat} = 1.0 \times M1$")
    fit3, = wp_ax.plot(wpsimfibre[:, 0], wpsimfibre[:, 1],
                       color='r', label=r"$M_{sat} = 1.5 \times M1$")
    fit4, = wp_ax.plot(wpsimfibre2[:, 0], wpsimfibre2[:, 1],
                       color='y', label=r"$M_{sat} = 2.0 \times M1$")

#    clust, = wp_ax.plot(clustering[:, 0], clustering[:, 1],
#                        color='y', label="analytical")

    wp_ax.errorbar(wpfile[:, 0], wpfile[:, 1], yerr=errors, fmt='.')

    wp_ax.set_xscale('log')
    wp_ax.set_yscale('log')

    handles, labels = wp_ax.get_legend_handles_labels()
    wp_ax.legend(handles, labels)


# plot residuals between fit and data

    fitspl = spline(wpsim[:, 0], wpsim[:, 1])

    res = wpfile[:, 1] - fitspl(wpfile[:, 0])
    res = res / wpfile[:, 1]

    reserrors = errors / wpfile[:, 1]

    res_ax.errorbar(wpfile[:, 0], res, reserrors)


# annotate

    wp_ax.set_title(r"$w_{p}$ fit")
    wp_ax.set_ylabel(r"$w_{p}(r_{p})$ ($h^{-1} Mpc$)")
#    wp_ax.text(0.6, 0.9, "{0} {1} {2}".format(*date), transform=wp_ax.transAxes)

    res_ax.set_xlabel(r"$r_{p}$ ($h^{-1} Mpc$)")
    res_ax.set_ylabel(r"Residuals")

    plotaddress = picname

    plt.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))


if __name__ == "__main__":

    if len(sys.argv) != 5:
        print "usage: python wp_errplot.py <wp datafile> <wp covar file> \
              <wp fit file> <png output filename>"

    plot_vpfs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
            sys.argv[6])

