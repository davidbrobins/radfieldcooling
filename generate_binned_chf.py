#Run as python generate_binned_chf.py <save_dir> <a_exp> <halo_mass> <num_halos>
# produces saved dataframes containing the binned cooling and heating functions
# save_dir (str): the path to the directory where the resulting dataframes should be saved
# a_exp (str): the scale factor a, as a string
# halo_mass (str): the desired mass of halos (in units of M_sun/h), should be a power of 10 (but doesn't need to be) OR the string 'most_massive'
# num_halos (str): the number of halos to bin (both individually, and aggregated), an integer

#Import commands 
import sys #Package to unpack inputs in call 
import os  #Package to confirm filepath exists
import numpy as np #Numpy (math)
import pandas as pd #Dataframe handling
from halo_catalog import HaloCatalog #Module defining HaloCatalog class
from halo_cooling_heating_functions import * #Package for cooling/heating function and cooling time binning

(pyfilename, save_dir, a_exp, halo_mass, num_halos)=sys.argv #Unpack the input variables (note: the element sys.argv[0], pyfilename, is the name of THIS FILE)

#Check that there is a directory in the path 'filepath'. If not, create it.

if not os.path.exists(save_dir+'/a'+a_exp+'/m'+halo_mass):
    os.makedirs(save_dir+'/a'+a_exp+'/m'+halo_mass+'/isochronous')
    os.makedirs(save_dir+'/a'+a_exp+'/m'+halo_mass+'/flowline')

#Read in the needed halo catalog (determined by the scale factor expression a_exp) and mass cuts
hc=HaloCatalog(a_exp) #Import the h5 halo catalog at the right scale factor
if halo_mass=='most_massive': #For 'most_massive', just take the entire dataframe
    halos_df=hc.halos_df
    maximum_mass=halos_df.iloc[0]['Mvir'] #Mass of highest-mass halo in the catalog
else: #If halo_mass is an integer, cut to virial masses between halo_mass, 1.2*halo_mass 
    halos_in_mass_bin=hc.select_halos(property='Mvir', min=float(halo_mass), max=1.2*float(halo_mass)) #Select the halos in the mass bin: between halo_mass and 1.2*halo_mass (might need tweaking)
    halos_df=halos_in_mass_bin.sort_values("Mvir", ascending=True) #Sort in ascending order (will take the first num_halos items
    maximum_mass=halos_df.iloc[int(num_halos)-1]['Mvir'] #Mass of highest-mass halo selected in this bin
#Write the highest virial mass, in units of Msun/h, to a file in the appropriate directory
max_mass_file=open(save_dir+'/a'+a_exp+'/m'+halo_mass+'/max_virial_mass.txt', 'w')
max_mass_file.write("Highest virial mass [Msun/h]: "+str(maximum_mass))
max_mass_file.close()

#Useful objects
field_list=['temperature', 'baryon_number_density', 'metallicity', ('artio', 'RT_DISK_VAR_0'), ('artio', 'RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio', 'RT_DISK_VAR_3'), 'cooling_rate', 'heating_rate', 'cooling_time'] #Needed fields for this analysis
all_halo_spheres_df=pd.DataFrame(columns=field_list) #Blank dataframe to hold field data for all 20 halos (concatenated)
all_halo_spheres_list=[] #List to hold all the dataframes (for looping)

#Loop through and get cell-by-cell data for each halo

for halo_index in range(int(num_halos)):
    df = get_halo_sphere_data(halos_df, "/nfs/turbo/lsa-cavestru/shared_data/CROC/dbrobins/rei20c1_a"+a_exp+"_RF/rei20c1_a"+a_exp+"_RF.art", halo_index, float(a_exp), field_list) #Import data
    df_cut=cut_fields_data(df, property='baryon_number_density', min=1, max=10) #Select cells with 1<n_b<10 [cm^-3]
    all_halo_spheres_list.append(df_cut) #Append the data frames to the list for later looping
    all_halo_spheres_df=all_halo_spheres_df.append(df_cut, ignore_index=True) #Concatenate all the halo spheres

#Set up temperature bins
all_halo_spheres_df=all_halo_spheres_df.sort_values(by='temperature') #Sort values by temperature (ascending)
bins=log_spaced_bins(all_halo_spheres_df['temperature'].to_numpy()[0], all_halo_spheres_df['temperature'].to_numpy()[-1], 20) #Create 20 log-spaced temperature bins between min and max temps

#Loop through halos and bin
for halo_index in range(int(num_halos)):
    df=all_halo_spheres_list[halo_Index] #Pick out the right halo sphere dataframe
    flowline_averaging(df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+halo_mass+'/flowline/halo_index_'+str(halo_index)) #Flowline binning
    isochronous_averaging(df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+halo_mass+'/isochronous/halo_index_'+str(halo_index)) #Isochronous binning

#Bin outside the halo
flowline_averaging(df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+halo_mass+'/flowline/stacked_halos') #Flowline binning
isochronous_averaging(df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+halo_mass+'/isochronous/stacked_halos') #Isochronous binning


