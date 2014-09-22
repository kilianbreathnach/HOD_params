import numpy as np
from astropy.cosmology import Planck13, WMAP5


def get_comv(cosmology):

    if cosmology == "Planck":
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif cosmology == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5

    return cosmo.comoving_distance


def get_inv_efunc(cosmology):

    if cosmology == "Planck":
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif cosmology == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5

    return cosmo.inv_efunc


