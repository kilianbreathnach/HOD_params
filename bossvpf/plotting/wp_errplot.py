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


if len(sys.argv) != 5:
    print "usage: python wp_errplot.py <wp datafile> <wp covar file> \
            <wp fit file> <png output filename>"

# Get date to add to plot
process = subprocess.Popen(['date', '+%d-%b-%Y'], stdout=subprocess.PIPE)
date, err = process.communicate()
date = date.split()[0].split("-")


# Get file data
wpfile = np.loadtxt(sys.argv[1])
clustering = np.loadtxt(sys.argv[3], usecols=(0, 4))

# Set up axis dimensions
wp_ax = plt.axes([0.1, 0.35, 0.9, 0.65])
res_ax = plt.axes([0.1, 0.1, 0.9, 0.2])


# plot wp from wp_covar from real mock catalogue against .clustering values

errors = np.log10(np.diag(np.loadtxt(sys.argv[2])))

wp_ax.errorbar(wpfile[:, 0], wpfile[:, 1], errors,
               linestyle='none', marker='o')
clust, = wp_ax.plot(clustering[:, 0], clustering[:, 1])

wp_ax.set_xscale('log')
wp_ax.set_yscale('log')


# plot residuals between fit and data

fitspl = spline(clustering[:, 0], clustering[:, 1])

res = wpfile[:, 1] - fitspl(wpfile[:, 0])
res = res / wpfile[:, 1]

reserrors = errors / wpfile[:, 1]

res_ax.errorbar(wpfile[:, 0], res, reserrors)
res_ax.set_xscale('log')


# annotate

wp_ax.set_title(r"$w_{p}$ fit")
wp_ax.set_ylabel(r"$w_{p}(r_{p})$ ($h^{-1} Mpc$)")
wp_ax.text(0.6, 0.9, "{0} {1} {2}".format(*date), transform=wp_ax.transAxes)

res_ax.set_xlabel(r"$r_{p}$ ($h^{-1} Mpc$)")
res_ax.set_ylabel(r"Residuals")

plt.savefig("/home/kilian/public_html/tinker/{0}.png".format(sys.argv[4]))
