import os
import numpy as np
from nbar_utils import nzobj


def process_sim(datrnk, simfile, sample, cap, cosmology, cutoffval):

    simout = datrnk + "/out/sim/" + sample + "/" + cap + "/" + cosmology +
    "/" + simfile

    sim_nz = nzobj(datrnk, sample, cap, cosmology, cutoffval)
    sim_nz.read_json()

    sim_nz.mock = np.genfromtxt(simout)

    sim_nz.fibre_collisions()

    hodroot = simfile.split(".")[0]

    print os.path.dirname(simout)
    sim_nz.sim_down(os.path.dirname(simout) + hodroot + ".down")
