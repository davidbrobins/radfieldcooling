import yt
from yt.analysis_modules.halo_analysis.api import *
import numpy as np
import pandas as pd
import h5py

h=0.6814000010490417
a=0.13557 
ds = yt.load("/data/dbrobins/40/C/rei40c1_a0.1356/rei40c1_a0.1356.art")
#df=pd.read_csv('/data/dbrobins/40/C/rs/hlists/hlist_0.13557.list', skiprows=64, header=None, delim_whitespace=True, names=['a', 'Mvir', 'Rvir', 'x', 'y', 'z', 'vx', 'vy', 'vz'], usecols=[0,10, 11, 17, 18, 19, 20, 21, 22], index_col=False) #Read in the halo from hlists
#most_massive_pandas = (df.sort_values('Mvir', ascending=False)).iloc[:10]
#print(most_massive_pandas)
#halos_ds = yt.load("/data/dbrobins/40/C/rs/halos_9.0.bin") #Have yt make a halo catalog from the binary rock star data
#hc = HaloCatalog(data_ds=ds, halos_ds=halos_ds)
#hc.create()
hl =yt.load('/home/dbrobins/repos/radfieldcooling/scripts/halo_catalogs/catalog/catalog.0.h5')
halos_bin = hl.all_data()
#print(halos_bin["halos","particle_mass"].in_units("Msun/h"))
df_bin = pd.DataFrame({"ID":halos_bin["halos", "particle_identifier"], "Mvir":halos_bin["halos","particle_mass"].in_units("Msun/h"), 'Rvir':halos_bin["halos", "virial_radius"].in_units("kpccm/h"), "x":halos_bin["halos","particle_position_x"].in_units("Mpccm/h"),"y":halos_bin["halos","particle_position_y"].in_units("Mpccm/h"),"z":halos_bin["halos","particle_position_z"].in_units("Mpccm/h")})
df_bin = df_bin.sort_values('Mvir', ascending=False) #Sort halos by virial mass (in descending order)
for i in range(3): #Pick the 3 most massive halos
    ID=df_bin.iloc[i]["ID"] #Halo ID#
    Mvir=df_bin.iloc[i]['Mvir'] #Halo virial mass
    R=df_bin.iloc[i]['Rvir'] #Halo virial radius
    x=df_bin.iloc[i]['x'] #Halo position
    y=df_bin.iloc[i]['y']
    z=df_bin.iloc[i]['z']
    halo_sphere = ds.sphere(yt.YTArray([x*a/h, y*a/h, z*a/h], "Mpc"),(2*R*a/h, "kpc"))
    halo_proj_density = yt.ProjectionPlot(ds, 1, 'density',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)                                                  
    halo_proj_temp = yt.ProjectionPlot(ds, 1, 'temperature',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)                                                 
    halo_proj_pressure = yt.ProjectionPlot(ds, 1, 'pressure',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)                                                
    halo_proj_entropy = yt.ProjectionPlot(ds, 1, 'entropy',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
    density_skewer = yt.LinePlot(ds, 'density', yt.YTArray([(x-2*R/1000)*a/h, y*a/h, z*a/h], "Mpc"), yt.YTArray([(x+2*R/1000)*a/h, y*a/h, z*a/h], "Mpc"), 512)
    halo_proj_density.save('/home/dbrobins/repos/radfieldcooling/plots/halo_%i_proj_density.png' % ID)                                                                                                
    halo_proj_temp.save('/home/dbrobins/repos/radfieldcooling/plots/halo_%i_proj_temp.png' % ID)
    halo_proj_pressure.save('/home/dbrobins/repos/radfieldcooling/plots/halo_%i_proj_pressure.png' % ID)                                                                                         
    halo_proj_entropy.save('/home/dbrobins/repos/radfieldcooling/plots/halo_%i_proj_entropy.png' % ID)            
    density_skewer.save('/home/dbrobins/repos/radfieldcooling/plots/halo_%i_skewer_density.png' % ID)                                                                                                    
#slice_of_halo = yt.SlicePlot(ds, 'z', 'density',center=([x*a/h,y*a/h,z*a/h], "Mpc"),
#                             width=(4*R*a/h, "kpc"))
#halo_proj_density = yt.ProjectionPlot(ds, 1, 'density',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
#halo_proj_temp = yt.ProjectionPlot(ds, 1, 'temperature',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
#halo_proj_pressure = yt.ProjectionPlot(ds, 1, 'pressure',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
#halo_proj_entropy = yt.ProjectionPlot(ds, 1, 'entropy',center = ([x*a/h, y*a/h, z*a/h], "Mpc"), width=(4*R*a/h, "kpc"), data_source=halo_sphere)
#slice_of_halo.set_width((2500, "kpccm/h"))
#halo_proj_density.save('/home/dbrobins/repos/radfieldcooling/plots/halo_proj_density_test.png')
#halo_proj_temp.save('/home/dbrobins/repos/radfieldcooling/plots/halo_proj_temp_test.png')
#halo_proj_pressure.save('/home/dbrobins/repos/radfieldcooling/plots/halo_proj_pressure_test.png')
#halo_proj_entropy.save('/home/dbrobins/repos/radfieldcooling/plots/halo_proj_entropy_test.png')
#slice_of_halo.save('/home/dbrobins/repos/radfieldcooling/plots/halo_slice_test.pdf')
#slice_of_halo.save('/home/dbrobins/repos/radfieldcooling/plots/halo_slice_test.png')
#density_skewer = yt.LinePlot(ds, 'density', yt.YTArray([(x-2*R/1000)*a/h, y*a/h, z*a/h], "Mpc"), yt.YTArray([(x+2*R/1000)*a/h, y*a/h, z*a/h], "Mpc"), 512)
#density_skewer.save('/home/dbrobins/repos/radfieldcooling/plots/skewer_density_test.png')

 
