from glob import glob
import numpy as np
from gen_plot import plotz

sph_files = glob("../1mil_hgg/volfrac*")

sph_count = []

for f in sph_files:
    sph_count.append(np.loadtxt(f))


plotz(masslist, deltarrs, style="hist_quant", axdim=[1, 3], xlog=True,
      xbins=len(sph_files), ybins=10, ymin=0.0, ymax=1.0
      name="/home/kilian/public_html/tinker/quantile_mass-delta.png")

