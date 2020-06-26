import yt
from yt.analysis_modules.halo_analysis.api import *
import numpy as np
import pandas as pd
import h5py

h=0.6814000010490417
a=0.1451 
ds = yt.load("/data/dbrobins/20/A/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art")
#df=pd.read_csv('/data/dbrobins/40/C/rs/hlists/hlist_0.13557.list', skiprows=64, header=None, delim_whitespace=True, names=['a', 'Mvir', 'Rvir', 'x', 'y', 'z', 'vx', 'vy', 'vz'], usecols=[0,10, 11, 17, 18, 19, 20, 21, 22], index_col=False) #Read in the halo from hlists
#most_massive_pandas = (df.sort_values('Mvir', ascending=False)
#halos_ds = yt.load("/data/dbrobins/40/C/rs/halos_9.0.bin") #Have yt make a halo catalog from the binary rock star data
#hc = HaloCatalog(data_ds=ds, halos_ds=halos_ds)
#hc.create()
#hl =yt.load('/home/dbrobins/repos/radfieldcooling/scripts/halo_catalogs/catalog/catalog.0.h5')
#halos_bin = hl.all_data()
#print(halos_bin["halos","particle_mass"].in_units("Msun/h"))
#df_bin = pd.DataFrame({"ID":halos_bin["halos", "particle_identifier"], "Mvir":halos_bin["halos","particle_mass"].in_units("Msun/h"), 'Rvir':halos_bin["halos", "virial_radius"].in_units("kpccm/h"), "x":halos_bin["halos","particle_position_x"].in_units("Mpccm/h"),"y":halos_bin["halos","particle_position_y"].in_units("Mpccm/h"),"z":halos_bin["halos","particle_position_z"].in_units("Mpccm/h")})
#df_bin = df_bin.sort_values('Mvir', ascending=False) #Sort halos by virial mass (in descending order)
x=0.825796
y=4.257736
z=4.266610
R=29.932201
halo_sphere = ds.sphere(yt.YTArray([x*a/h, y*a/h, z*a/h], "Mpc"),(2*R*a/h, "kpc"))
print(dir(ds.fields.artio))
P_LW= yt.SlicePlot(ds, 'z',('artio', 'RT_DISK_VAR_0'),center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)                                                  
P_LW.set_colorbar_label('RT_DISK_VAR_0', '$P_{LW}$ (1/s)')
P_HI= yt.SlicePlot(ds, 'z', ('artio','RT_DISK_VAR_1'),center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
P_HI.set_colorbar_label('RT_DISK_VAR_1', '$P_{HI}$ (1/s)')
P_HeI= yt.SlicePlot(ds, 'z',('artio', 'RT_DISK_VAR_2'),center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
P_HeI.set_colorbar_label('RT_DISK_VAR_2', '$P_{HeI}$ (1/s)')
P_CVI= yt.SlicePlot(ds, 'z',('artio', 'RT_DISK_VAR_3'),center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
P_CVI.set_colorbar_label('RT_DISK_VAR_3', '$P_{CVI}$ (1/s)')
#P_CVI.set_zlim(field='RT_DISK_VAR_3',zmin=10e-19, zmax=10e-17)
dens= yt.SlicePlot(ds, 'z', 'density',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
temp= yt.SlicePlot(ds, 'z', 'temperature',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
P_LW.save('/home/dbrobins/repos/radfieldcooling/plots/rad_field_tests_4_6/P_LW.png')
P_HI.save('/home/dbrobins/repos/radfieldcooling/plots/rad_field_tests_4_6/P_HI.png')                                                                                         
P_HeI.save('/home/dbrobins/repos/radfieldcooling/plots/rad_field_tests_4_6/P_HeI.png')
P_CVI.save('/home/dbrobins/repos/radfieldcooling/plots/rad_field_tests_4_6/P_CVI.png')
dens.save('/home/dbrobins/repos/radfieldcooling/plots/rad_field_tests_4_6/dens.png')
temp.save('/home/dbrobins/repos/radfieldcooling/plots/rad_field_tests_4_6/temp.png')
