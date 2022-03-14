#Run as python generate_binned_chf.py <save_dir> <a_exp> <halo_mass> <num_halos>
# produces saved dataframes containing the binned cooling and heating functions
# save_dir (str): the path to the directory where the resulting dataframes should be saved
# a_exp (str): the scale factor a, as a string
# min_halo_mass (str): the desired minimum mass of halos (in units of M_sun/h) [DO NOT GO BELOW 1e10 for resolution reasons]
# num_to_bin (str): The number of halos to randomly subselect and bin (if available) 
# halo_catalog_dir (str): the directory path to the .list halo catalogs
# central_or_sub (str): 'central' or 'sub' to denote central or subhalos, respectively
# ind_or_all (str): 'ind' or 'all' to denote whether halos should also be binned one-by-one ('ind') or just combined prior to binning ('all')


#Import commands 
import sys #Package to unpack inputs in call 
import os  #Package to confirm filepath exists
import numpy as np #Numpy (math)
import pandas as pd #Dataframe handling
from hlist_halo_catalog import hlistHaloCatalog #Module defining HaloCatalog class
from halo_cooling_heating_functions import * #Package for cooling/heating function and cooling time binning

(pyfilename, save_dir, a_exp, min_halo_mass, num_to_bin, halo_catalog_dir, central_or_sub, ind_or_all)=sys.argv #Unpack the input variables (note: the element sys.argv[0], pyfilename, is the name of THIS FILE)

#Check that there is a directory in the path 'filepath'. If not, create it.

if not os.path.exists(save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/'):
    os.makedirs(save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/isochronous')
    os.makedirs(save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/flowline')

#Read in the needed halo catalog (determined by the scale factor expression a_exp) and mass cuts
hc=hlistHaloCatalog(a_exp, halo_catalog_dir=halo_catalog_dir) #Read in the hlist halo catalog at the right scale factor
halos_df=hc.halos_df #Get the dataframe
if central_or_sub=='central': #If want central halos, cut to just central halos
    halos_df=halos_df.loc[halos_df['pid']==-1]
elif central_or_sub=='sub': #If want subhalos, cut to just subhalos
    halos_df=halos_df.loc[halos_df['pid']>-1]
halos_in_mass_bin= halos_df.loc[(halos_df['Mvir']>float(min_halo_mass))] #Select the halos above the minimum mass
halos_df=halos_in_mass_bin.sort_values("Mvir", ascending=True) #Sort in ascending order (will take the first num_halos items)
num_halos=str(len(halos_df.index))
print(num_halos)

#Useful objects
field_list=['temperature', 'baryon_number_density', 'metallicity', ('artio', 'RT_DISK_VAR_0'), ('artio', 'RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio', 'RT_DISK_VAR_3'), 'cooling_rate', 'heating_rate', 'cooling_time'] #Needed fields for this analysis
all_halo_spheres_df=pd.DataFrame(columns=field_list) #Blank dataframe to hold field data for all 20 halos (concatenated)
all_halo_spheres_list=[] #List to hold all the dataframes (for looping)
    
total_cells=0 # counter to check the total number of (gas) cells after the density cut
kept_cells=0 # counter to check how many cells survive the radiation field cuts
    
# Check if there are enough halos to randomly subselect
if int(num_to_bin) < int(num_halos): # if so
    halo_indices = np.random.choice(int(num_halos), int(num_to_bin), replace = False) # randomly choose num_to_bin halo indices
else: # otherwise
    halo_indices = range(int(num_halos)) # choose all halos


# Loop through the list of halo indices defined above
for halo_index in halo_indices: 
    df = get_halo_sphere_data(halos_df, "/nfs/turbo/lsa-cavestru/shared_data/CROC/dbrobins/rei20c1_a"+a_exp+"_RF/rei20c1_a"+a_exp+"_RF.art", halo_index, float(a_exp), field_list) #Import data for the halo
    df_cut=cut_fields_data(df, property='baryon_number_density', min=1, max=10) #Select cells with 1<n_b<10 [cm^-3]
    total_cells += len(df_cut.index) # Add the number of cells to total cell number
    # Cut cells with P_{HI, HeI, CVI} <= 0 (CHF approximation doesn't work there)
    final_df=df_cut.loc[(df_cut[('artio','RT_DISK_VAR_1')] > 0.0) & (df_cut[('artio','RT_DISK_VAR_2')] > 0.0) & (df_cut[('artio','RT_DISK_VAR_3')] > 0.0)]
    kept_cells += len(final_df.index) # Add the number of surviving cells to surviving cell number
    all_halo_spheres_list.append(final_df) #Append the data frames to the list for later looping
    all_halo_spheres_df=all_halo_spheres_df.append(final_df, ignore_index=True) #Concatenate all the halo spheres

#Set up temperature bins
all_halo_spheres_df=all_halo_spheres_df.sort_values(by='temperature') #Sort values by temperature (ascending)
bins=log_spaced_bins(all_halo_spheres_df['temperature'].to_numpy()[0], all_halo_spheres_df['temperature'].to_numpy()[-1], 20) #Create 20 log-spaced temperature bins between min and max temps
    
discarded_fraction = (total_cells-kept_cells) / total_cells

# Write the halo indices used to a .txt file for future reference
info_file = open(save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/'+str(len(halo_indices))+central_or_sub+'_halos_info.txt', 'w')
info_file.write('Number of halos: '+str(len(halo_indices))+'\n')
info_file.write('Halo indices: \n')
for halo_index in halo_indices:
    info_file.write('%i \n' % halo_index)
info_file.write('Cell fraction discarded in radiation field cut: %.3f' % discarded_fraction)
info_file.close()

# Do the binning
flowline_averaging(all_halo_spheres_df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/flowline/'+central_or_sub+'_halos') #Flowline binning 
print('Did flowline binning')
isochronous_averaging(all_halo_spheres_df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/isochronous/'+central_or_sub+'_halos') #Instantaneous binning 
print('Did instantaneous binning')
if ind_or_all == 'ind': # If 'ind', bin the halos individually
    for index in range(len(halo_indices)): # Now just loop through number of halos that we're using
        halo_df=all_halo_spheres_list[index] #Pick out the right halo sphere dataframe
        flowline_averaging(halo_df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/flowline/'+central_or_sub+'_'+str(index)) #Flowline binning # remove the number of halos later                                                                                           
        print('Did flowline binning for halo '+str(index)+' of '+str(len(halo_indices)))
        isochronous_averaging(halo_df, bins, binning_field='temperature', to_save=True, filepath=save_dir+'/a'+a_exp+'/m'+min_halo_mass+'/isochronous/'+central_or_sub+'_'+str(index)) #Instantaneous binning # remove the number of halos later                                                                               
        print('Did instantaneous binning for halo '+str(index)+' of '+str(len(halo_indices)))
