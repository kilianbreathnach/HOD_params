from glob import glob
from dat.h5_funcs import mk_h5, fits2h5
from dat.nbar_zrange import mk_zfile


sdss_dir = "./dat/in/sdss3"
CMASS_dir = "./dat/in/CMASS_DATA"
out_dir = "./dat/out/CMASS/NGC"
cosmo = "WMAP"
cosmo_dir = "{0}/{1}".format(out_dir, cosmo)

mk_zfile("{0}/nbar-cmass-dr11v1-N-Anderson.dat".format(sdss_dir),
         "{0}/full_radecz.hdf5".format(out_dir), cosmo)


