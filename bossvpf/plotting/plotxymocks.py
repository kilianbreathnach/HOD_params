import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

mock = np.loadtxt("../../dat/out/sim/dzwide.mock", usecols=(0, 1))
mockdown = np.loadtxt("../../dat/out/sim/simdown.dat", usecols=(0, 1))

ax1.scatter(mock[:, 0], mock[:, 1],
                  alpha=0.05)
ax2.scatter(mockdown[:, 0], mockdown[:, 1],
                 alpha=0.05)

picname = "plotxy.png"
fig.savefig(picname)

os.system("scp {0} broiler:~/public_html/tinker/".format(picname))
os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))



