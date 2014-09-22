import json
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, Distance, ICRSCoordinates
from astropy.cosmology import Planck13, WMAP5
from scipy.spatial import cKDTree
from dat.h5_funcs import h5_arr

"""
This script takes a file of galaxy positions that have been transformed from
BOSS data to a cartesian coordinate system using a given cosmology. It sorts
them into python's k-d tree and then drops sphere's randomly in the space to
find the probability for void spheres at various radii.

The random points are constructed from a file of random angular coordinates
distributed in the BOSS mask according to completeness and from a redshift
taken from a uniform sample in a spherical shell. This redshift shell was
determined from the number density of observed galaxies in the survey, taken
to be within an interval of the maximum of the number density distribution.

The spheres are then tested for their integrity by checking how much of their
volume lies outside of the mask. This is done by using a sample of points with
given solid angle density. These are also sorted into a k-d tree of their
locations on the unit sphere. By finding the number of points lying within
the circular projection of the  void sphere onto the unit 2-sphere, and based
on their distance from the centre of the circle, the fractional volume of the
sphere whose voidness is uncertain may be determined.
"""


# Define function to help build k-d tree for bad points and do search
def radec2xyz(radecarr):

    radecarr = np.atleast_2d(radecarr)
    xyzarr = np.zeros((radecarr.shape[0], 3))
    xyzarr[:, 0] = np.cos(np.radians(radecarr[:, 1])) * \
                   np.cos(np.radians(radecarr[:, 0]))
    xyzarr[:, 1] = np.cos(np.radians(radecarr[:, 1])) * \
                   np.sin(np.radians(radecarr[:, 0]))
    xyzarr[:, 2] = np.sin(np.radians(radecarr[:, 1]))

    return xyzarr

# And a function to find angular distance of bad points to sphere
# projection centres
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


def vpf(dat_dir, Nsph, simul_cosmo, rad):

    # Grab the data coordinates
    gals = h5_arr("./dat/out/{0}/{1}/gals_cart_coords.hdf5".
                      format(dat_dir, simul_cosmo), "cart_pts")

    # Get details about the redshift interval being considered
    nbar_dict = json.load(open("./dat/out/{0}/{1}/nbar_zrange.json".
                                   format(dat_dir, simul_cosmo)))
    zlo = nbar_dict["zlo"]
    zhi = nbar_dict["zhi"]

    # Get the search points
    good_pts = h5_arr("./dat/out/{0}/srch_radec.hdf5".format(dat_dir), "good_pts")
    bad_pts = h5_arr("./dat/out/{0}/veto.hdf5".format(dat_dir),
                     "bad_pts")

    # Set angular radius of effective area around bad points
    bad_r = np.arccos(1.0 - (np.pi * 9.8544099e-05) / (2 * 180 ** 2))
    bad_r_deg = np.rad2deg(bad_r)

    # Set the cosmology with h free
    # Here the cosmology is based on WMAP (for first MultiDark simulation)
    if simul_cosmo == "Planck":
        # First make h free
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif simul_cosmo == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5
    comv = cosmo.comoving_distance

    # Build the trees

    # galaxy tree
    gal_baum = cKDTree(gals)

    # tree of bad points (angular coordinates on unit sphere)
    bad_xyz = radec2xyz(bad_pts)
    veto_baum = cKDTree(bad_xyz)

    # Initialise final output arrays
#    rad = np.arange(1.0, 67.0, 5.0)  doing it one radius at a time
#    P_0 = np.zeros(rad.shape)

    # No. of spheres and norm
#     Nsph_arr = Nsph * np.array(4 * [0.01] + 4 * [0.1] + 4 * [1.0])
#     norm = 1. / Nsph_arr
#    norm = 1. / Nsph

    rand_i = 0

    for r_i, r in enumerate(rad):

        # start the count of successful voids
        count = 0

        # Custom zrange for sphere size
        dis_near = Distance(comv(zlo).value + r, u.Mpc)
        dis_far = Distance(comv(zhi).value - r, u.Mpc)

        z_a = dis_near.compute_z(cosmology=cosmo)

        z_b = dis_far.compute_z(cosmology=cosmo)

        for i in range(Nsph):  # _arr[r_i]):

            # compensate for finite length of mask file
            rand_i = rand_i % 999999

            radec = good_pts[rand_i, :]

            rang = Angle(radec[0], u.deg)
            decang = Angle(radec[1], u.deg)

            randz = (z_a ** 3 + \
                     (z_b ** 3 - z_a ** 3) * np.random.rand(1)[0]) ** (1. / 3.)
            dis = Distance(comv(randz), u.Mpc)

            coord = ICRSCoordinates(rang, decang, distance=dis)

            sph_cen = np.array([coord.x.value, coord.y.value, coord.z.value])

            nn = gal_baum.query(sph_cen)

            print "rad: ", r, ", sphere: ", i

            if not nn[0] < r:

                # add instance to probability count
                count += 1

                # record quality of sphere using spline values for intersection
                # with bad points

                # Get radius of circular projection of sphere
                R = np.arcsin(r / np.sqrt(np.sum(sph_cen[:] ** 2)))

                # Get coordinates of circle centre on unit sphere
                crc_cen = radec2xyz(radec)[0]

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
                    dis = np.degrees(central_angle(pt_ang, radec))
                    l = dis / R

                    bad_vol += 1.5 * (bad_r_deg / R) ** 2 \
                                   * np.sqrt(1.0 - l ** 2)

                f_r = open("./dat/out/{0}/{1}/vpf_out/volfrac_{2}.dat".
                               format(dat_dir, simul_cosmo, r),
                           'a')
                f_r.write("{0}\n".format(bad_vol))
                f_r.close()

            rand_i += 1

        # Save the void probability regardless of integrity of void
        # (can decide on cutoff for reliability afterwards)
#        P_0[r_i] = norm * count  # [r_i] * count

    # Print VPF to file
#    fin_out = np.zeros((rad.shape[0], 2))

#    fin_out[:, 0] = rad
#    fin_out[:, 1] = P_0

#    np.savetxt("./dat/{0}/{1}/out/rad-P_0.dat".format(dat_dir, simul_cosmo), fin_out)


if __name__ == "__main__":

    import sys

    # Get the runtime arguments

    if len(sys.argv) != 4:
        print "usage: python vpf_analysis.py \
               <survey-cap dirname: [CMASS, LOWZ]/[NGC, SGC] OR 'theoretical'>\
               <No. of spheres at each radius> <Planck | WMAP>"

    vpf(sys.argv[1], int(sys.argv[2]), sys.argv[3])
