import json
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck13, WMAP5
from dat.h5_funcs import h5_arr


Nfrac = (6769.0358 * np.pi) / 129600

WMAP5.__init__(100.0, WMAP5.Om0)
cosmo = WMAP5
comv = cosmo.comoving_distance

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


nz_f = open("dat/out/CMASS/NGC/WMAP/nbar_zrange.json")
nbar_file = "dat/in/sdss3/nbar-cmass-dr11v1-N-Anderson.dat"
nz_dict = json.load(nz_f)
nbar_corr = np.loadtxt(nbar_file)
nbar_corr = nbar_corr[(nz_dict["zlo"] <= nbar_corr[:, 0]) * \
                      (nbar_corr[:, 0] <= nz_dict["zhi"])]
zcen = nbar_corr[:, 0]
z_near = nbar_corr[:, 1]
z_far = nbar_corr[:, 2]

print "zcen: ", zcen
print "z_near: ", z_near
print "z_far: ", z_far

dzn = z_near[1] - z_far[0]
dzf = z_far[1] - z_far[0]

z_near = np.arange(z_near[0], z_near[-1] + dzn, 0.5 * dzn)
z_far = np.arange((0.5 * (z_near[0] + z_near[1])), z_far[-1] + 0.5 * dzf, 0.5 * dzf)
zcen = 0.5 * (z_far + z_near)
zbinedges = np.append(z_near[0], z_far[:])

print "\nthen\n"
print "zcen: ", zcen
print "z_near: ", z_near
print "z_far: ", z_far
print "zbinedges: ", zbinedges

shell_vols = []
for i in range(len(zcen)):
    shell_vols.append(Nfrac * calc_shell_vol(comv, z_near[i], z_far[i], zcen[i]))
shell_vols = np.array(shell_vols)

dat_rdz = h5_arr("dat/out/CMASS/NGC/WMAP/radecz_down.hdf5", "radecz")

dat_hist = np.histogram(dat_rdz[:, 2], bins=zbinedges)[0]

dat_nbar = dat_hist / shell_vols

mock_nbars = []
mock_fns = glob("dat/out/CMASS/NGC/WMAP/mocks/rdz_down/*")

for mock in mock_fns:

    radecz = h5_arr(mock, "radecz")

    H = np.histogram(radecz[:, 2], bins=zbinedges)

    mock_nbars.append(H[0] / shell_vols)

# and plotting

fig = plt.figure()
ax = fig.add_subplot(111)

datl, = ax.plot(zcen, dat_nbar, color="r")

mockls = []
for mock in mock_nbars:

    mockls.append(ax.plot(zcen, mock, color = "b", alpha=0.5))

fig.savefig("downsamples.pdf")


