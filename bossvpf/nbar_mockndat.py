import numpy as np
from astropy.io import fits


datf = fits.open("CMASS/NGC/cmass-dr11v1-N-Anderson.dat.fits")

nbar_arr = np.loadtxt("CMASS/NGC/nbar-cmass-dr11v1-N-Anderson.dat")

mock_zs = np.loadtxt("/data1/kww231/mocks/ngc/a0.6452_0292.dr11c_ngc.rdz", usecols=(2,))
dat_zs = datf[1].data["Z"]

zbinedges = np.append(nbar_arr[0, 1], nbar_arr[:, 2])

zcens = nbar_arr[:, 0]

dat_corr = nbar_arr[:, 6]
mockH = np.histogram(mock_zs, bins=zbinedges)
datH = np.histogram(dat_zs, bins=zbinedges)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure()

ax = fig.add_subplot(111)

l_corr, = ax.plot(zcens, dat_corr, label="corrected number density")
l_mock, = ax.plot(zcens, mockH[0], label="mock number density")
l_dat, = ax.plot(zcens, datH[0], label="real data")

handles, labels = ax.get_legend_handles_labels()

ax.legend(handles, labels)

fig.savefig("compare_nbar_mock.pdf")
