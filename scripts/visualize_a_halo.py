import yt
import numpy as np
import pandas as pd

ds = yt.load("/data/dbrobins/C/rei20c1_a0.1352/rei20c1_a0.1352.art")
halo_sphere = ds.sphere([7.22082, 15.82621, 2.38917],(7500, "kpc/h"))

slice_of_halo = yt.SlicePlot(ds, 'z', 'density', center=[7.22082, 15.82621, 2.38917], width=(7500, "kpc/h"))
slice_of_halo.set_width((7500, "kpc/h"))
slice_of_halo.save('/home/dbrobins/repos/radfieldcooling/plots/halo_slice_test.pdf')
slice_of_halo.save('/home/dbrobins/repos/radfieldcooling/plots/halo_slice_test.png')

