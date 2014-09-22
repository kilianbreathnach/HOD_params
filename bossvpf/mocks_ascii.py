import os
import json
import numpy as np
from scipy.spatial import cKDTree
from astropy import units as u
from astropy.coordinates import Angle, Distance, ICRSCoordinates
from astropy.cosmology import Planck13, WMAP5
from vpf_calc import vpf
from dat.cartesian_cosmo import mk_coords
from dat.h5_funcs import h5_arr, arr2h5


def mk_mock_srch(radecfile, nzdictfile, Nsph, simul_cosmo):

    if simul_cosmo == "Planck":
        # First make h free
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif simul_cosmo == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5
    comv = cosmo.comoving_distance

    radecarr = h5_arr(radecfile, "good_pts")
    nzdict = json.load(open(nzdictfile))

    Nrands = radecarr.shape[0]
    Narrs = Nsph / Nrands
    remain = Nsph % Nrands

    radecz = np.zeros((Nsph, 3))

    for i in range(Narrs):

        start = Nrands * i
        stop = Nrands * (i + 1)
        radecz[start:stop, :2] = radecarr[:, :]

    endchunk = Nrands * (Narrs)
    radecz[endchunk:, :2] = radecarr[:remain, :]

    rad = np.arange(1.0, 67.0, 5.0)
    zlo = nzdict["zlo"]
    zhi = nzdict["zhi"]

    radeczlist = len(rad) * [radecz]

    for r_i, r in enumerate(rad):

        dis_near = Distance(comv(zlo).value + r, u.Mpc)
        dis_far = Distance(comv(zhi).value - r, u.Mpc)

        z_a = dis_near.compute_z(cosmology=cosmo)

        z_b = dis_far.compute_z(cosmology=cosmo)

        randz = (z_a ** 3 + \
                 (z_b ** 3 - z_a ** 3) * np.random.rand(Nsph)) ** (1. / 3.)

        radeczlist[r_i][:, 2] = randz[:]

        arr2h5(radeczlist[r_i], "{0}/{1}/mocks/mock_srch_pts.hdf5".format(os.path.dirname(radecfile), simul_cosmo), "radecz_{0}".format(str(r_i * 5 + 1)))


def mk_mock_coords(radeczfile, outfile, simul_cosmo):

    if simul_cosmo == "Planck":
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif simul_cosmo == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5

    rad = np.arange(1.0, 67.0, 5.0)

    radecz = h5_arr(radeczfile, "radecz")

    cart = np.zeros(radecz.shape)

    for i, rdz in enumerate(radecz):

        ra = Angle(rdz[0], u.deg)
        dec = Angle(rdz[1], u.deg)

        losd = cosmo.comoving_distance(rdz[2])
        dis = Distance(losd)

        coord = ICRSCoordinates(ra, dec, distance=dis)

        cart[i, :] = np.array([coord.x.value, coord.y.value, coord.z.value])

    np.savetxt(outfile, cart)


def mock_vpf(mock_cart_coords, spheresfile, simul_cosmo, rad):

    gals = h5_arr(mock_cart_coords, "coords")

    print gals

    name = mock_cart_coords.split("/")[-1].split(".")[0]

    if simul_cosmo == "Planck":
        # First make h free
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif simul_cosmo == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5
    comv = cosmo.comoving_distance

    gal_baum = cKDTree(gals)

    spheres = h5_arr(spheresfile, "radecz_{0}".format(str(int(rad))))

    print spheres

    for i, sphere in enumerate(spheres):

        rang = Angle(sphere[0], u.deg)
        decang = Angle(sphere[1], u.deg)

        dis = Distance(comv(sphere[2]), u.Mpc)

        coord = ICRSCoordinates(rang, decang, distance=dis)

        sph_cen = np.array([coord.x, coord.y, coord.z])

        nn = gal_baum.query(sph_cen)

        print "rad: ", rad, ", sphere: ", i

        f = open("{0}/vpf_out/ascii/{1}_{2}.dat".format(os.path.dirname(spheresfile), name, str(int(rad))), 'a')

        if not nn[0] < rad:
            f.write("1\n")
        else:
            f.write("0\n")

        f.close()
