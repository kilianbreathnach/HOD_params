import os
import json
import random
from glob import glob
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from cosmo_stuff import get_comv
from scipy.spatial import cKDTree


def calc_shell_vol(dis_func, z_near, z_far, zcen):
    """
    Computes the volume of a spherical shell of at redshift zcen, with edges
    at redshifts z_near and z_far, given a cosmological distance conversion
    function. Results are in (h^-1 Mpc)^3.
    """

    r = dis_func(zcen).value
    r_near = dis_func(z_near).value
    r_far = dis_func(z_far).value
    dr = r_far - r_near

    return 4 * np.pi * (r ** 2) * dr


def polynom(x, *coeffs):
    res = np.zeros(x.shape)
    for i in range(len(coeffs)):
        res += coeffs[i] * x ** i
    return res


class nzobj:

    def __init__(self, dattop, sample, cap, cosmology, cutoff):

        self.dattop = dattop
        self.sample = sample
        self.cap = cap

        if cap == "NGC":
            self.shellfrac = (6769.0358 * np.pi) / 129600

        self.cosmology = cosmology
        self.cutoff = cutoff

        self.comv = get_comv(cosmology)

        self.nzjson = dattop + "/out/{0}/{1}/{2}/nbar_zrange_{3}.json".format(sample, cap, cosmology, cutoff)


    def get_zrange(self):

        # corrected number density file for processing
        nbarfile = glob(self.dattop + "/in/boss/{0}/{1}/nbar*.dat".format(self.sample,
            self.cap))[0]
        self.nbar_corr = np.genfromtxt(nbarfile)

        # Cut out the first bit of crap (works for CMASS, dunno about LOWZ)
        ind03 = np.abs(self.nbar_corr[:, 0] - 0.3).argmin()
        self.nbar_corr = self.nbar_corr[ind03:, :]

        zcen = self.nbar_corr[:, 0]
        z_near = self.nbar_corr[:, 1]
        z_far = self.nbar_corr[:, 2]
        corr_gal_counts = self.nbar_corr[:, 6]

        nbar = []
        shell_vols = []

        for i in range(len(zcen)):

            shell_vols.append(self.shellfrac * \
                              calc_shell_vol(self.comv, z_near[i], z_far[i], zcen[i]))
            nbar.append(corr_gal_counts[i] / shell_vols[i])

        nbar = np.array(nbar)
        shell_vols = np.array(shell_vols)

        # Find nbar peak and index
        max_nbar = np.max(nbar)
        max_i = int(np.where(nbar == max_nbar)[0])

        self.nz_dict = {"tophat height for zrange": self.cutoff}
        self.nz_dict["max_nbar_corr"] = max_nbar
        self.nz_dict["nbar_corr_tophat"] = self.cutoff * max_nbar
        self.nz_dict["z_nbar_max"] = zcen[max_i]

        # get the interval edge indices
        L = np.abs(nbar[:max_i] - max_nbar * self.cutoff).argmin()
        R = max_i + np.abs(nbar[max_i:] - max_nbar * self.cutoff).argmin()

        self.nbar_full = nbar
        self.nbar = nbar[L:R + 1]
        self.zcen = zcen[L:R + 1]
        shell_vols = shell_vols[L:R + 1]

        self.nz_dict["zlo"] = zcen[L]
        self.nz_dict["zhi"] = zcen[R]
        self.nz_dict["dz"] = zcen[1] - zcen[0]

        self.nz_dict["avg_nbar_corr"] = np.average(nbar)
        self.nz_dict["total_shell_vol"] = np.sum(shell_vols)


    def read_json(self):

        jf = open(self.nzjson)
        self.nz_dict = json.load(jf)
        jf.close()


    def _histy(self, buf=0):

        z0 = self.nz_dict["zlo"] - (buf + .5) * self.nz_dict["dz"]
        zn = self.nz_dict["zhi"] + (buf + .5) * self.nz_dict["dz"]

        self.z_near = np.arange(z0, zn, self.nz_dict["dz"])
        self.z_cen = self.z_near + self.nz_dict["dz"] / 2.
        self.z_far = self.z_near + self.nz_dict["dz"]

        self.zbinedges = np.append(z0, self.z_far)
        self.hist = np.histogram(self.radecz[:, 2], bins=self.zbinedges)[0]


    def _get_nbar(self):

        nbar = []
        shell_vols = []

        for i in range(len(self.z_cen)):
            shell_vols.append(self.shellfrac * \
                calc_shell_vol(self.comv, self.z_near[i],
                               self.z_far[i], self.z_cen[i]))

        shell_vols = np.array(shell_vols)
        self.nbar = self.hist / shell_vols



    def _compute_z_nbar(self, buf=0, histonly=False):

        self._histy(buf)

        if not histonly:
            self._get_nbar()

            if buf == 0:
                self.nz_dict["avg_nbar_dat"] = np.average(self.nbar)


    def fit_polynom(self):

        # we're fitting 4th order polynomials here
        deg = 4

        self._compute_z_nbar(buf=5)
        self.popt, self.pcov = curve_fit(polynom, self.z_cen, self.nbar,
                                         p0=np.random.rand(deg))
        result = polynom(self.z_cen, *self.popt)

        # save the results
        np.savetxt(os.path.dirname(self.nzjson) + "/poly_fit.dat",
                   np.dstack((self.z_cen, result))[0])
        np.savetxt(os.path.dirname(self.nzjson) + "/poly_coeffs.dat", self.popt)


    def rdz_down(self, outfile):

        self._compute_z_nbar()

        polyfit = polynom(self.z_cen, *self.popt)

        # check if we're dealing with a mock here
        ismock = False
        outsplit = outfile.split("/")
        if outsplit[1] == "mocks":
            ismock = True

        print "ismock:\t", ismock

        # mocks are to be sampled down to min from data
        if ismock:
            ndown = self.nz_dict["data_tophat"]
        else:
            ndown = np.min(self.nbar)
            self.nz_dict["data_tophat"] = ndown

        print "ndown:\t", ndown

        # downsampling ratio
        keepfrac = ndown / self.nbar

        # for mocks, we shouldn't have a keepfrac greater than one
        keepfrac[np.where(keepfrac > 1.0)[0]] = 1

        print "keepfrac"
        print keepfrac

        # Get the raw histogram
        self._compute_z_nbar(histonly=True)

        # Get an integer number to downsample to
        num_down = np.rint(keepfrac * self.hist)
        num_down = num_down.astype(int)

        print "numdown"
        print num_down

        finmask = np.array(self.radecz.shape[0] * [True])

        for i, nd in enumerate(num_down):
            """Turn on the right amount of galaxies in each bin."""
            zbin_ids = np.where(((self.zbinedges[i] < self.radecz[:, 2]) *\
                (self.radecz[:, 2] <= self.zbinedges[i + 1])) == True)

            bincount = zbin_ids[0].shape[0]
            if bincount == 0:
                continue

            Ngo = bincount - nd

            if Ngo <= 0:
                print i, Ngo
                continue

            toss = np.random.choice(zbin_ids[0], size=Ngo, replace=False)

            finmask[toss] = False

        self.radecz = self.radecz[finmask]

        np.savetxt(os.path.dirname(self.nzjson) + outfile,
                   self.radecz)


    def _xydist(xy1, xy2):
        return np.sqrt((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2)


    def fibre_collisions(self):

        # redshift of simulation snapshot
        z = 0.5342

        # Fibre collision completeness
        thres = 0.42

        # Transverse comoving distance of fibre at simulation redshift
        fibrad = np.radians(62. / 3600) * self.comv(z).value

        # take out the xy coords
        mockxy = self.mock[:, :2]

        # set up a tree :)
        gal_baum = cKDTree(mockxy)

        # Search for points inside fibre radius of each point
        collisions = gal_baum.query_ball_tree(gal_baum, fibrad)

        colmask = np.array(mockxy.shape[0] * [True])

        # randomise how the galaxies are checked for collisions
        indlist = range(len(collisions))
        random.shuffle(indlist)

        coldict = {}

        for ind in indlist:

            if colmask[ind] == False:
                continue

            if len(collisions[ind]) > 1:

                coldict[ind] = [self.mock[a].tolist() for a in
                                collisions[ind]]

#                 for galind in collisions[ind]:
#                     if np.random.rand() > thres:
#                         colmask[galind] = False

                # don't kill the galaxy iterated over
                collisions[ind].pop(collisions[ind].index(ind))

                for galind in collisions[ind]:
                    colmask[galind] = False

                # randomly choose one of the collided galaxies and
                # keep it if it beats probability threshold
                if np.random.rand() <= thres:
                    galind = random.choice(collisions[ind])
                    colmask[galind] = True


        fcol = open("../dat/out/sim/collisions.dat", 'w')
        json.dump(coldict, fcol, sort_keys=True, indent=4, separators=(',', ':\t'))

        fcol.close()

        origlen = len(self.mock)

        self.mock = self.mock[colmask]

        newlen = len(self.mock)

        print "downsampled by collisions to ", float(newlen) / origlen


    def sim_down(self, outfile):

        box_nbar = float(len(self.mock)) / (1000 ** 3)
        print box_nbar
        keepfrac = self.nz_dict["data_tophat"] / box_nbar
        print self.nz_dict["data_tophat"]

        # number of galaxies to sample down to
        numdown = np.rint(keepfrac * len(self.mock))
        numdown = numdown.astype(int)
        N2go = len(self.mock) - numdown

        finmask = np.array(self.mock.shape[0] * [True])
        mockinds = np.array(range(len(self.mock)))

        toss = np.random.choice(mockinds, size=N2go, replace=False)

        finmask[toss] = False

        self.mock = self.mock[finmask]

        print "final sim nbar: ", float(len(self.mock)) / (1000 ** 3)

        np.savetxt(outfile, self.mock)


    def write_json(self):

        nf = open(self.nzjson, 'w')

        json.dump(self.nz_dict, nf, sort_keys=True, indent=4, separators=(',', ':\t'))

        nf.close()
