import os
import sys
from glob import glob
import numpy as np



def mock_frac_vars(dirnam):

    allmocks = glob(dirnam + "/*.vpf")

    mockvpfs = []

    for fn in allmocks:

        mockvpfs.append(np.genfromtxt(fn, usecols=(1,)))

    vpfarr = []

    for mock in mockvpfs:
        try:
            if not mock.shape[0] < 8:
                vpfarr.append(mock)
        except:
            pass

    mockvpfs = np.array(vpfarr)

    zs = np.loadtxt(allmocks[0], usecols=(0,))

    meanP = []

    for i in range(len(zs)):

        meanP.append(np.mean(mockvpfs[:, i]))

    meanP = np.array(meanP)


    diffsq = np.zeros(mockvpfs.shape)

    for i in range(len(zs)):
        diffsq[:, i] = (mockvpfs[:, i] - meanP[i]) ** 2

    meandiffsq = np.zeros(meanP.shape)

    for i in range(len(zs)):
        meandiffsq[i] = np.mean(diffsq[:, i])

    print meandiffsq
    frac_var = meandiffsq / (meanP ** 2)
    print frac_var

    np.savetxt(dirnam + "/frac_var.dat", np.dstack((zs, meanP, meandiffsq, frac_var))[0])



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "plot_vpf takes exactly 1 argument"

    mock_frac_vars(sys.argv[1])
