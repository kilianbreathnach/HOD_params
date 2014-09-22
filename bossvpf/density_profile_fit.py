import os
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


nz_dict = {}
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

# Save values
nz_dict["max_nbar_corr"] = max_nbar
nz_dict["nbar_corr_tophat"] = Q * max_nbar
nz_dict["z_nbar_max"] = zcens[max_i]

# get the interval edge indices
L = np.abs(nbar[:max_i] - max_nbar * Q).argmin()
R = max_i + np.abs(nbar[max_i:] - max_nbar * Q).argmin()

fit_nbar = nbar[L - 5:R + 6]
nbar = nbar[L:R + 1]
shell_vols = shell_vols[L:R + 1]

nz_dict["zlo"] = zcens[L]
nz_dict["zhi"] = zcens[R]

nz_dict["avg_nbar_corr"] = np.average(nbar)
nz_dict["total_shell_vol"] = np.sum(shell_vols)



# Do this shizz for the data with different bin widths!

# first with the given widths
z0 = z_near[L - 5]
zn = z_far[R + 5]

nbar = get_n_vols(dat_zs, z0, zn, dz)

# then for smaller bin widths
nbar2 = get_n_vols(dat_zs, z0, zn, dz * 0.1)


# Now compare these in a plot!

gs = grd.GridSpec(4, 4, wspace=0.8, hspace=0.8)

# Set up axis dimensions
n_ax = plt.subplot(gs[:3, :2])
n2_ax = plt.subplot(gs[:3, 2:])

res_ax = plt.subplot(gs[-1, :2])
res2_ax = plt.subplot(gs[-1, 2:])

z_near = np.arange(z0, zn, dz)
z_cen = z_near + dz / 2.

z2_near = np.arange(z0, zn, dz * 0.1)
z2_cen = z2_near + (dz * 0.1) / 2.

# get plotting
n_ax.scatter(z_cen, nbar, alpha=0.5, label="nbar values")
n2_ax.scatter(z2_cen, nbar2, alpha=0.5, label="nbar values")

for degree in [3, 4, 5]:
    popt, pcov = curve_fit(polynom, z_cen, nbar, p0=np.random.rand(degree))
    result = polynom(z_cen, *popt)
    n_ax.plot(z_cen, result, alpha=0.5, label="degree %d" % degree)
    resid = nbar - result
    res_ax.plot(z_cen, resid, alpha=0.5, label="degree %d" % degree)

    popt, pcov = curve_fit(polynom, z2_cen, nbar2, p0=np.random.rand(degree))
    result2 = polynom(z2_cen, *popt)
    n2_ax.plot(z2_cen, result2, alpha=0.5, label="degree %d" % degree)
    resid2 = nbar2 - result2
    res2_ax.plot(z2_cen, resid2, alpha=0.5, label="degree %d" % degree)


n_ax.set_ylim(0, 0.0005)
n_ax.legend(loc=1, prop={'size':8})
res_ax.legend(loc=1, prop={'size':4})

n2_ax.set_ylim(0, 0.0005)
n2_ax.legend(loc=1, prop={'size':8})
res2_ax.legend(loc=1, prop={'size':4})

picname = "polynomial_nbar_fits.pdf"

plt.savefig(picname)

os.system("scp {0} broiler:~/public_html/tinker/".format(picname))
os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
              .format(picname))



"""
# Now compute histogram for the data

nbar_arr = nbar_arr[(nz_dict["zlo"] <= nbar_arr[:, 0]) * \
                     (nbar_arr[:, 0] <= nz_dict["zhi"])]

# number density histogram bin edges
zbinedges = np.append(nbar_arr[0, 1], nbar_arr[:, 2])

H = np.histogram(dat_zs, bins=zbinedges)
dat_nbar = H[0] / shell_vols




# Now to fit a polynomial to the corrected number densities

# # create matrix versions of these arrays
# zcens = zcens[L - 5:R + 6]
# X = zcens[:, np.newaxis]
#
# plt.scatter(zcens, fit_nbar, label="nbar values")
#
# for degree in [3, 4, 5]:
#     model = make_pipeline(PolynomialFeatures(degree), Ridge())
#     model.fit(X, fit_nbar)
#     y_plot = model.predict(X)
#     plt.plot(zcens, y_plot, alpha=0.5, label="degree %d" % degree)
#
# plt.legend(loc='lower left')
#
# plt.savefig("polynomial_nbar_fits.pdf")

zcens = zcens[L - 5:R + 6]


# Set up axis dimensions
n_ax = plt.axes([0.1, 0.35, 0.9, 0.65])
res_ax = plt.axes([0.1, 0.1, 0.9, 0.2])


n_ax.scatter(zcens, fit_nbar, alpha=0.5, label="nbar values")

# for degree in [3, 4, 5]:
#     p = np.polyfit(zcens, fit_nbar, degree, rcond=None, full=False, w=None, cov=False)
#     result = np.zeros(zcens.shape)
#     for i in range(degree):
#         result += p[i] * zcens ** i
#     plt.plot(zcens, result, alpha=0.5, label="degree %d" % degree)



# now try the scipy way

def polynom(x, *coeffs):
    res = np.zeros(x.shape)
    for i in range(len(coeffs)):
        res += coeffs[i] * x ** i
    return res

for degree in [3, 4, 5]:
    popt, pcov = curve_fit(polynom, zcens, fit_nbar, p0=np.random.rand(degree))
    result = polynom(zcens, *popt)
    n_ax.plot(zcens, result, alpha=0.5, label="degree %d" % degree)
    resid = fit_nbar - result
    res_ax.plot(zcens, resid, alpha=0.5, label="degree %d" % degree)


n_ax.set_ylim(0, 0.0005)
n_ax.legend(loc=1)
res_ax.legend(loc=1, prop={'size':6})

picname = "polynomial_nbar_fits.pdf"

plt.savefig(picname)

os.system("scp {0} broiler:~/public_html/tinker/".format(picname))
os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
              .format(picname))

"""
