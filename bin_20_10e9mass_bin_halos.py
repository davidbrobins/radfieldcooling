#A script to get isochronous and flowline cooling and heating functions for 20 halos in the 10^9 Msun/h mass bin, binned separately and together
import numpy as np #numpy
import pandas as pd #pandas (dataframe handling)
from halo_catalog import HaloCatalog #HaloCatalog() class
from halo_cooling_heating_functions import *

hc=HaloCatalog("/home/dbrobins/repos/radfieldcooling/radfieldcooling/scripts/halo_catalogs/RF_0.1451_catalog/RF_0.1451_catalog.0.h5") #open the needed halo catalog 
halos_df=hc.halos_df #Extract the dataframe of the halo properties
halos_mass_10e9_df=halos_df.loc[(halos_df['Mvir'] >= 10**9) & (halos_df['Mvir'] <= 1.1*10**9)] #Cut just the halos in [1, 1.1]*10^9 Msun/h mass bin
all_50_halo_spheres_df=pd.DataFrame(columns=['temperature', 'baryon_number_density','metallicity', ('artio', 'RT_DISK_VAR_0'), ('artio', 'RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio', 'RT_DISK_VAR_3'), 'cooling_rate', 'heating_rate', 'cooling_time']) #Blank dataframe to hold field data for all 20 halos (concatenated)
list_of_halo_spheres=[] #Blank list to hold the each dataframe
for halo_index in range(50): #Loop through 20 halos
    df = get_halo_sphere_data(halos_mass_10e9_df, "/nfs/turbo/lsa-cavestru/shared_data/CROC/dbrobins/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art", halo_index, 0.1451, ['temperature', 'baryon_number_density','metallicity', ('artio', 'RT_DISK_VAR_0'), ('artio', 'RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio', 'RT_DISK_VAR_3'), 'cooling_rate', 'heating_rate', 'cooling_time'])
    df_cut = cut_fields_data(df, property='baryon_number_density', min=1, max=10)                                                         
    list_of_halo_spheres.append(df_cut) #append the dataframes to a list for later looping
    all_50_halo_spheres_df=all_50_halo_spheres_df.append(df_cut, ignore_index=True)
all_50_halo_spheres_df=all_50_halo_spheres_df.sort_values(by='temperature')
bins=log_spaced_bins(all_50_halo_spheres_df['temperature'].to_numpy()[0], all_50_halo_spheres_df['temperature'].to_numpy()[-1], 20) #Create 20 log spaced bins per dex in T between min. and max. temps for all cells in all 50 halos (with 1<n_b<10 cm^-3)
flowline_averaging(all_50_halo_spheres_df, bins, binning_field='temperature', to_save=True, name='50_at_mass_10e9')
isochronous_averaging(all_50_halo_spheres_df, bins, binning_field='temperature', to_save=True, name='50_at_mass_10e9')
for halo_index in range(50):
    df=list_of_halo_spheres[halo_index]
    flowline_df = flowline_averaging(df, bins, binning_field='temperature', to_save=True, name='mass_10e9_'+str(halo_index))
    isochronous_df=isochronous_averaging(df, bins, binning_field='temperature', to_save=True, name='mass_10e9_'+str(halo_index))
