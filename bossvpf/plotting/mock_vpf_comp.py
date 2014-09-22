from glob import glob
import matplotlib.pyplot as plt
import numpy as np


mockfiles = glob("./dat/out/tinkercode/output/mocks/????.out")

datfile = glob("./dat/out/tinkercode/output/vpf/10percent.dat")
dat_arr = np.loadtxt(datfile[0])
zcen = dat_arr[:, 0]
dat_vpf = np.log10(dat_arr[:, 1])


# plot

fig = plt.figure()
ax = fig.add_subplot(111)

mockls = []
for mock in mockfiles:

    try:
        mock_vpf = np.log10(np.loadtxt(mock)[:, 1])
        mockls.append(ax.plot(zcen, mock_vpf, color = "b", alpha=0.05))
    except:
        continue

datl, = ax.plot(zcen, dat_vpf, color="r")

fig.savefig("vpf_mock_compare.pdf")


