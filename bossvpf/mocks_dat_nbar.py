import os
from glob import glob
import numpy as np
import matplotlib
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
from scipy.optimize import curve_fit

from astropy.io import fits
from astropy.coordinates import Distance
from astropy.cosmology import Planck13, WMAP5

from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline


def calc_shell_vol(dis_func, z_near, z_far, zcen):
    """
    Computes the volume of a spherical shell of at redshift zcen, with edges
    at redshifts z_near and z_far, given a cosmological distance conversion
    function. Results are in (h^-1 Mpc)^3.
    """

    r = dis_func(zcen).value
    r_near = dis_func(z_near).value
    r_far = dis_func(z_far).value
    dr = r_far - r_near

    return 4 * np.pi * (r ** 2) * dr


def get_n_vols(datzs, z0, zn, dz):

    z_near = np.arange(z0, zn, dz)
    z_cen = z_near + dz / 2.
    z_far = z_near + dz

    shell_vols = []
    nbar = []

    for i in range(len(z_cen)):

        shell_vols.append(shellfrac * calc_shell_vol(comv, z_near[i],
                                                     z_far[i], z_cen[i]))

    shell_vols = np.array(shell_vols)
    zbinedges = np.append(z_near[0], z_far)
    H = np.histogram(datzs, bins=zbinedges)
    nbar = H[0] / shell_vols

    return nbar


def polynom(x, *coeffs):
    res = np.zeros(x.shape)
    for i in range(len(coeffs)):
        res += coeffs[i] * x ** i
    return res


# magic number for width around maximum
Q = 0.65
# magic number for shell vol computation (in NGC)
shellfrac = (6769.0358 * np.pi) / 129600

cosmology = "WMAP"

if cosmology == "Planck":
    Planck13.__init__(100.0, Planck13.Om0)
    cosmo = Planck13
elif cosmology == "WMAP":
    WMAP5.__init__(100.0, WMAP5.Om0)
    cosmo = WMAP5
comv = cosmo.comoving_distance


# load the corrected density data into memory
nbar_arr = np.loadtxt("../dat/in/boss/CMASS/NGC/nbar-cmass-dr11v1-N-Anderson.dat")

# Cut out the first bit of crap (works for CMASS, dunno about LOWZ)
ind03 = np.abs(nbar_arr[:, 0] - 0.3).argmin()
nbar_arr = nbar_arr[ind03:, :]

# get the galaxy data to compute number density from
galfits = fits.open("../dat/in/boss/CMASS/NGC/cmass-dr11v1-N-Anderson.dat.fits")
# array of the redshifts
dat_zs = galfits[1].data["Z"]

# redshift arrays (bin edges and centres)
zcens = nbar_arr[:, 0]
z_near = nbar_arr[:, 1]
z_far = nbar_arr[:, 2]

dz = 0.005

# corrected number densities at each redshift
corr_counts = nbar_arr[:, 6]

# get mock files and their redshifts
mockfns = glob("../dat/in/boss/CMASS/NGC/mocks/a0.6452_00[0-2]*.rdz")

mockzs = []
for fn in mockfns:
    mockzs.append(np.loadtxt(fn, usecols=(2,)))

nbar = []
shell_vols = []

for i in range(len(zcens)):

    shell_vols.append(shellfrac * calc_shell_vol(comv, z_near[i], z_far[i], zcens[i]))
    nbar.append(corr_counts[i] / shell_vols[i])

nbar = np.array(nbar)
shell_vols = np.array(shell_vols)

# Find nbar peak and index
max_nbar = np.max(nbar)
max_i = int(np.where(nbar == max_nbar)[0])

# get the interval edge indices
L = np.abs(nbar[:max_i] - max_nbar * Q).argmin()
R = max_i + np.abs(nbar[max_i:] - max_nbar * Q).argmin()

fit_nbar = nbar[L - 5:R + 6]
nbar_corr = nbar[L:R + 1]


# Do this shizz for the data and mocks with different bin widths!

# first with the given widths
z0 = z_near[L - 5]
zn = z_far[R + 5]

nbar = get_n_vols(dat_zs, z0, zn, dz)

factor = nbar_corr / nbar[5:-5]
bigF = np.mean(factor)

nbar = bigF * nbar

mockbars = []
for zs in mockzs:
    mockbars.append(get_n_vols(zs, z0, zn, dz))

# then for smaller bin widths
nbar2 = bigF * get_n_vols(dat_zs, z0, zn, dz * 0.1)

mock2bars = []
for zs in mockzs:
    mock2bars.append(get_n_vols(zs, z0, zn, dz * 0.1))


# Now compare these in a plot!

gs = grd.GridSpec(4, 4, wspace=0.8)

# Set up axis dimensions
n_ax = plt.subplot(gs[:, :2])
n2_ax = plt.subplot(gs[:, 2:])

z_near = np.arange(z0, zn, dz)
z_cen = z_near + dz / 2.

z2_near = np.arange(z0, zn, dz * 0.1)
z2_cen = z2_near + (dz * 0.1) / 2.

# get plotting
n_ax.scatter(z_cen, nbar, color='r', alpha=0.5, label="data nbar values")
for mock in mockbars:
    n_ax.plot(z_cen, mock, color='b', alpha=0.3)

n2_ax.scatter(z2_cen, nbar2, color='r', alpha=0.5, label="data nbar values")
for mock in mock2bars:
    n2_ax.plot(z2_cen, mock, color='b', alpha=0.3)

n_ax.set_ylim(0, 0.0005)
n_ax.legend(loc=1, prop={'size':8})

n2_ax.set_ylim(0, 0.0005)
n2_ax.legend(loc=1, prop={'size':8})

picname = "mock_dat_nbar.pdf"

plt.savefig(picname)

os.system("scp {0} broiler:~/public_html/tinker/".format(picname))
os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
              .format(picname))
