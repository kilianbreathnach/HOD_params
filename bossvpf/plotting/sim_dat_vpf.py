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

#     golden_mean = (np.sqrt(5)-1.0)/2.0
#     w, h = plt.figaspect(float(3 / golden_mean))
#
#     fig = plt.figure(figsize=(w,h))

#     fig, axes = plt.subplots(ncols=3, figsize=(10, 6),
#                              subplot_kw=dict(aspect=1, adjustable='box'))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    colours = ['g', 'm', 'r', 'y']
    lstyles = [":", "--", "-.", '-']
    lines = []

    vpf = np.loadtxt(dattop + "/out/sim/subhalos.vpf")
    lines.append(ax.plot(np.append(0., vpf[:, 0]),
                         np.append(1., vpf[:, 1]),
                         linestyle=lstyles.pop(), color=colours.pop(),
                         alpha=0.5,
                         label=r"subhalos mock"))

    vpf2 = np.loadtxt(dattop + "/out/CMASS/NGC/WMAP/vpf/sim_comp_0.95.vpf")
    lines.append(ax.plot(np.append(0., vpf2[:, 0]),
                         np.append(1., vpf2[:, 1]),
                         linestyle=lstyles.pop(), color=colours.pop(),
                         alpha=0.5,
                         label=r"normal mock"))

#     for i, frac in enumerate(["0.5", "1.0", "1.5", "2.0"]):
#
#         simfn = dattop + "/out/CMASS/NGC/WMAP/vpf/Mtest/{0}_42.vpf".format(frac)
#
#         vpf = np.loadtxt(simfn)
#         lines.append(ax.plot(np.append(0., vpf[:, 0]),
#                              np.append(1., vpf[:, 1]),
#                              linestyle=lstyles.pop(), color=colours.pop(),
#                              alpha=0.5,
#                              label=r"$M_{{sat}} = {0} \times M1$".format(frac)))
#
#         sim2fn = dattop + "/out/CMASS/NGC/WMAP/vpf/Mtest/{0}_42.vpf".format(frac)
#
#         vpf = np.loadtxt(sim2fn)
#         l2, = ax.plot(np.append(0., vpf[:, 0]),
#                       np.append(1., vpf[:, 1]),
#                       linestyle='--', color='m', alpha=0.5, label="42\%")
#
#         sim3fn = dattop + "/out/CMASS/NGC/WMAP/vpf/Mtest/{0}_down.vpf".format(frac)
#
#         vpf3 = np.loadtxt(sim3fn)
#         l2, = ax.plot(np.append(0., vpf3[:, 0]),
#                       np.append(1., vpf3[:, 1]),
#                       linestyle='-.', color='r', alpha=0.5, label="100\%")

    datfn = dattop + "/out/CMASS/NGC/WMAP/vpf/comp_0.05.vpf"
    vpfd = np.loadtxt(datfn)

#    l2, = ax.plot(np.append(0., vpf[:, 0]),
#                  np.append(1., vpf[:, 1]),
#                  marker='o', color='b',
#                  alpha=0.5, label='data')

    errors = np.loadtxt(dattop + "/out/CMASS/NGC/WMAP/vpf/mocks/frac_var.dat", usecols=(1,))

    ax.errorbar(np.append(0., vpfd[:, 0]),
                np.append(1., vpfd[:, 1]),
                yerr=np.append(0., vpfd[:, 1] * np.sqrt(errors)),
                fmt='.', lw=1, capsize=5, mew=1)
#                linestyle='none', marker='.', fmt='')
#                 marker='o', color='k', ecolor='k',
#                 markerfacecolor='b', label='data')
     #            linestyle=None) # fmt=None, color='b')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc=3, fontsize=10)

    ax.set_yscale('log', nonposy='clip')

    ax.set_xlabel(r"$r$ $(h^{-1} Mpc)$")
    ax.set_ylabel(r"$P_{0}(r)$")
#    ax.set_title(r"$M_{{sat}} = {0} \times (1.145051 \times 10^{{14}})$".format(frac))


    picname = "subhalos_vs_normal.pdf"
    plotaddress = os.path.dirname(datfn) + "/plots/" + picname

#    plt.tight_layout()

    fig.savefig(plotaddress)

    os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
    os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                  .format(picname))


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument: the address of the data directory"


    plot_vpfs(sys.argv[1])
