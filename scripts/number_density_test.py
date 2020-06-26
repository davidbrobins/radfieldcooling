import yt
from yt.analysis_modules.halo_analysis.api import *
from yt.units import gram
import numpy as np
import pandas as pd
import h5py

amu=1.66e-24*gram #atomic mass unit in g
def _mass_density_to_n_b(field,data): 
    return 0.26*data["density"]/(4.002602*amu)*4+(1-0.26)*data["density"]/(1.00784*amu)*1
h=0.6814000010490417
a=0.1451 
ds = yt.load("/data/dbrobins/20/A/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art")
ds.add_field(('gas', 'baryon_number_density'), function=_mass_density_to_n_b, units='1/cm**3')
x=0.825796
y=4.257736
z=4.266610
R=29.932201
halo_sphere = ds.sphere(yt.YTArray([x*a/h, y*a/h, z*a/h], "Mpc"),(2*R*a/h, "kpc"))
dens= yt.SlicePlot(ds, 'z', 'density',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
number_dens=yt.SlicePlot(ds, 'z', 'baryon_number_density', center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
dens.save('/home/dbrobins/repos/radfieldcooling/plots/density_comparison_tests_5_21/mass_density.png')
number_dens.save('/home/dbrobins/repos/radfieldcooling/plots/density_comparison_tests_5_21/baryon_number_density.png')
