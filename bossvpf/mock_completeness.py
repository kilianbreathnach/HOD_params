import os
import json
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, Distance, ICRSCoordinates
from astropy.cosmology import Planck13, WMAP5
from scipy.spatial import cKDTree
from h5_funcs import h5_arr, arr2h5

survey_cap = "CMASS/NGC"
simul_cosmo = "WMAP"

spheresfile = "out/{0}/{1}/mocks/mock_srch_pts.hdf5".format(survey_cap, simul_cosmo)

bad_pts = np.loadtxt("in/CMASS_DATA/north_block_outside.dat", usecols=(0, 1))

nbar_vals = json.load(open("out/{0}/{1}/nbar_zrange.json".format(survey_cap, simul_cosmo)))
zlo = nbar_vals["zlo"]
zhi = nbar_vals["zhi"]

bad_r = np.arccos(1.0 - (np.pi * 9.8544099e-05) / (2 * 180 ** 2))
bad_r_deg = np.rad2deg(bad_r)

if simul_cosmo == "Planck":
    Planck13.__init__(100.0, Planck13.Om0)
    cosmo = Planck13
elif simul_cosmo == "WMAP":
    WMAP5.__init__(100.0, WMAP5.Om0)
    cosmo = WMAP5
comv = cosmo.comoving_distance


def radec2xyz(radecarr):

    radecarr = np.atleast_2d(radecarr)
    xyzarr = np.zeros((radecarr.shape[0], 3))
    xyzarr[:, 0] = np.cos(np.radians(radecarr[:, 1])) * \
                   np.cos(np.radians(radecarr[:, 0]))
    xyzarr[:, 1] = np.cos(np.radians(radecarr[:, 1])) * \
                   np.sin(np.radians(radecarr[:, 0]))
    xyzarr[:, 2] = np.sin(np.radians(radecarr[:, 1]))

    return xyzarr


def central_angle(coord1, coord2):

    f1 = np.radians(coord1[1])
    f2 = np.radians(coord2[1])

    dl = abs(coord1[0] - coord2[0])

    if dl > 180:
        dl = np.radians(dl - 180)
    else:
        dl = np.radians(dl)
    # angle is (from wikipedia formula...)
    dsig = np.arccos(np.sin(f1) * np.sin(f2) + \
                     np.cos(f1) * np.cos(f2) * np.cos(dl))

    return dsig


bad_xyz = radec2xyz(bad_pts)
veto_baum = cKDTree(bad_xyz)

rad = np.arange(1.0, 67.0, 5.0)

rand_i = 0

for r_i, r in enumerate(rad):

    spheres = h5_arr(spheresfile, "radecz_{0}".format(str(r_i * 5 + 1)))

    badvols = np.zeros(spheres.shape[0])

    for i, sphere in enumerate(spheres):

        rang = Angle(sphere[0], u.deg)
        decang = Angle(sphere[1], u.deg)

        dis = Distance(comv(sphere[2]), u.Mpc)

        coord = ICRSCoordinates(rang, decang, distance=dis)

        sph_cen = np.array([coord.x.value, coord.y.value, coord.z.value])

        print "rad: ", r, ", sphere: ", i

        # Get radius of circular projection of sphere
        R = np.arcsin(r / np.sqrt(np.sum(sph_cen[:] ** 2)))

        # Get coordinates of circle centre on unit sphere
        crc_cen = radec2xyz(sphere[:2])[0]

        # Compute tree search radius from Cosine rule
        # (include points extending beyond sphere edge to account for
        # finite area around bad points)
        l_srch = np.sqrt(2. - 2. * np.cos(R))

        # Run search
        pierce_l = veto_baum.query_ball_point(crc_cen, l_srch)

        bad_vol = 0.

        R = np.degrees(R)  # need in degrees for bad_vol computation

        for pt in pierce_l:

            pt_ang = bad_pts[pt]
            dis = np.degrees(central_angle(pt_ang, sphere[:2]))
            l = dis / R

            bad_vol += 1.5 * (bad_r_deg / R) ** 2 * np.sqrt(1.0 - l ** 2)

        badvols[i] = bad_vol

    arr2h5(badvols,
            "{0}/mock_badvols.hdf5".format(os.path.dirname(spheresfile)),
            "badvols_{0}".format(str(r_i * 5 + 1)))
