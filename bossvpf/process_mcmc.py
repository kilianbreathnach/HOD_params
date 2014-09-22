import os
from glob import glob
import numpy as np
from nbar_utils import nzobj


def process_mcmc(datrnk, sample, cap, cosmology, cutoffval):


    mcmcfns = glob("{0}/out/sim/{1}/{2}/{3}/mcmc/new/*.mock".
                                 format(datrnk, sample, cap, cosmology))

    sim_nz = nzobj(datrnk, sample, cap, cosmology, cutoffval)
    sim_nz.read_json()

    downfns = glob("{0}/out/sim/{1}/{2}/{3}/mcmc/new/*.down".
                                 format(datrnk, sample, cap, cosmology))


    for mock in mcmcfns:

        simdown = mock[:-4] + "down"
        print simdown

        if not simdown in downfns:

            try:
                sim_nz.mock = np.genfromtxt(mock)
                sim_nz.fibre_collisions()
                sim_nz.sim_down(simdown)
            except:
                continue
