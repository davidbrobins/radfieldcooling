# A module to do binning and create isochronous and flowline averaged cooling/heating function arrays for a given halo
import yt #We need yt to do this
import numpy as np #Numpy
import pandas as pd #Pandas (dataframe handling)
import CF #f2py module that contains python wrappers for calculating cooling and heating functions
from halo_catalog import HaloCatalog #Module defining HaloCatalog() class
h=0.6814000010490417 #The hubble constant: H0=100*h km/s/Mpc
def log_spaced_bins(min, max, bins_per_dex): 
    '''
    Function to calculate flowline-averaged cooling and heating functions (and cooling time)
    Inputs: 
    min (float or int): The left side of the lowest bin
    max (float or int): The right side of the hightest bin
    bins_ber_dex (int): The number of bins per factor of 10 
    Outputs:
    bins (ndarray): An 1*(bins_per_dex*np.int(np.log10(max/min))) array of the bin edges
    '''
    num_bins=np.int(np.log10(max/min))*bins_per_dex #Find the number of log-spaced bins
    bins=np.logspace(np.log10(min), np.log10(max), num_bins) #Create the log-spaced bins
    return bins

def get_bin_centers(bins): 
    '''
    Function to get the centers (in log space) of log-spaced bins
    Inputs:
    bins (ndarray): a 1*n array of the bin edges (including left-most and right-most edges)
    Outputs:
    centers (ndarray): a 1*(n-1) array of the centers (in log space) of the bins
    '''
    centers = 10**((np.log10(bins[1:])+np.log10(bins[:-1]))/2) #Find the centers of each bin
    return centers

def get_halo_sphere_data(halos_df, ds_filename, halo_index, a, fields):
    '''
    Function to extract a dataframe with given fields and halo info dictionary for a given halo index
    Inputs:
    halos_df (dataframe): A dataframe with the given halos, likely imported from HaloCatalog() class
    ds (ytObject): The yt data set being used
    halo_index (int): The index of the halo in halos_df 
    a (float): scale factor [might be removed]
    fields (array of strings): Fields to be selected.  Choose from 'baryon_number_density', 'temperature', 'metallicity', ('artio', 'RT_DISK_VAR_[0-3]'), 'cooling_rate', 'heating_rate', 'cooling_time'
    Outputs:
    df (dataframe): Dataframe containing fields as the column labels
    '''
    halo_info=halos_df.iloc[halo_index] #Extract a dictionary with halo ID, Mvir, Rvir, x, y, z
    import derived_fields_ch_nb  #Import module with baryon_number_density, cooling_rate, heating_rate, cooling_time derived fields
    ds = yt.load(ds_filename)
    halo_sphere=ds.sphere(yt.YTArray([halo_info['x']*a/h, halo_info['y']*a/h, halo_info['z']*a/h], "Mpc"), (halo_info['Rvir']*a/h, "kpc")) #Pick out a Rvir-radius sphere centered on the halo center
    df=halo_sphere.to_dataframe(fields) #Extract a dataframe with the needed fields
    return df

def cut_fields_data(df, property='baryon_number_density', min=1, max=10):
    '''
    Function to get dataframe like df from get_halo_sphere_data to only between min-max values of property
    Inputs:
    df (dataframe): df like that from get_halo_sphere_data
    property (string): A valid field in df
    min: minimum value of property
    max: maximum value of property
    Outputs:
    cut_df (dataframe): df with rows where property<min or property>max removed
    '''
    cut_df=df.loc[(df[property]>min)&(df[property]<max)]
    return cut_df

def percentile_25(array):
    '''
    Function to return the 25th percentile of an array
    Inputs: 
    array: an array or ndarray
    Outputs:
    val (float): the 25th percentile of array
    '''
    return np.percentile(array, 25)

def percentile_75(array):
    #Function to return the 75th percentile of an array
    return np.percentile(array, 75)

def flowline_averaging(df, bins, binning_field='temperature'): #add an argument about saving later savenpz=True
    '''
    Function to take a dataframe and given bins for binning_field and return flowline averages (and 25th, 75th percentiles) for cooling_rate, heating_rate, cooling_time in those bins
    Inputs:
    df (dataframe): The dataframe used for averaging, must (at minimum) contain columns of binning_field, cooling_rate, heating_rate, and cooling_time
    bins (ndarray): Array of (logspaced) bins including left and right edges
    binning_field (string): The field to bin in (default to temperature)
    TO DO to_save (bool): Whether or not to save the resulting arrays (defaults to True)
    Outputs: 
    flowline_averages_df (dataframe): dataframe with centers of bin in binning_field and median, 25th, 75th percentile for cooling rate, heating rate, and cooling time
    '''
    centers=get_bin_centers(bins) #Get (log-spaced) centers of the bins
    flowline_averages_df = pd.DataFrame(columns=['bin_centers', 'CF_median', 'CF_25', 'CF_75', 'HF_median', 'HF_25', 'HF_75', 'ctime_median', 'ctime_25', 'ctime_75'])
    for bin_index in range(len(bins)-1): #Loop through the left edges of each bin
        binned_df=cut_fields_data(df, property=binning_field, min=bins[bin_index], max=bins[bin_index+1]) #Cut df to only contain the data where binning_field is within the appropriate bin
        cols_for_averaging=binned_df[['cooling_rate', 'heating_rate', 'cooling_time']] #Extract the three columns we want to average
        if len(cols_for_averaging.index) > 0: #Check the cols_for_averaging has at least one column.  If so, take the median and 25th, 75th percentiles of each column
            cf_med, hf_med, ct_med = cols_for_averaging.apply(np.median, 0) #Get the median of each column and store
            cf_25, hf_25, ct_25 = cols_for_averaging.apply(percentile_25, 0) #Get the 25th percentile of each column
            cf_75, hf_75, ct_75 = cols_for_averaging.apply(percentile_75, 0) #Get the 75th percentile of each column
            flowline_averages_df=flowline_averages_df.append({'bin_centers':centers[bin_index],'CF_median':cf_med, 'CF_25':cf_25, 'CF_75':cf_75, 'HF_median':hf_med, 'HF_25':hf_25, 'HF_75':hf_75, 'ctime_median':ct_med, 'ctime_25':ct_25, 'ctime_75':ct_75}, ignore_index=True) #Append the resulting values to the output dataframe
    return flowline_averages_df #Return the dataframe after the loop completes
    #if savenpz==True: save centers, c_func_median, h_func_median, c_time_median, percentiles

def CHF_from_array(input_array): #get cooling function from array of input values (for use in derived field calculation)                                                                         
    '''
    Function to get cooling and heating functions from an input array
    Input:
    input_array (array): Array containing [temperature, baryon number density, metallicity, P_LW (RT_DISK_VAR_0), P_HI (RT_DISK_VAR_1), P_HeI (RT_DISK_VAR_2), P_CVI (RT_DISK_VAR_3)]
    Outputs:
    cfun, hfun (floats): The cooling and heating functions at those values
    '''
    T=input_array[0]
    n_b=input_array[1]
    Z=input_array[2]
    P_LW=input_array[3]
    P_HI=input_array[4]
    P_HeI=input_array[5]
    P_CVI=input_array[6] #extract the 7 input floats from the array                                                                                                                                        
    (cfun,hfun,ierr)=CF.frtgetcf(T, n_b, Z, P_LW, P_HI, P_HeI, P_CVI) #call f2py function to get cooling and heating functions for these inputs                                                            
    return [cfun, hfun] #Return the cooling and heating functions

#def isochronous_averaging(df, bins, binning_field='temperature'): #add a boolean argument to_save that defaults to True
    #get centers of bins with centers=get_bin_centers(bins)
    #binning_field should be either T ('temperature') or n_b  ('baryon_number_density')
    #loop through centers, add a column to df which is just that T value
    #to df['T_for_calc', 'baryon_number_density', 'metallicity', 'RT_DISK_VAR_0', 'RT_DISK_VAR_1', 'RT_DISK_VAR_2', 'RT_DISK_VAR_3'].apply()
    #Take median of the result
    #if savenpz==True: save the results
    

#testing:

hc=HaloCatalog("/home/dbrobins/repos/radfieldcooling/scripts/halo_catalogs/RF_0.1451_catalog/RF_0.1451_catalog.0.h5")
halos_df=hc.halos_df
df = get_halo_sphere_data(halos_df, "/data/dbrobins/20/A/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art", 1, 0.1451, ['temperature', 'baryon_number_density','metallicity', ('artio', 'RT_DISK_VAR_0'), ('artio', 'RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio', 'RT_DISK_VAR_3'), 'cooling_rate', 'heating_rate', 'cooling_time'])
#print(df)
df_cut = cut_fields_data(df) 
#print(df_cut)
bins=log_spaced_bins(1, 1000, 20)
centers=get_bin_centers(bins)
#flowline_df = flowline_averaging(df_cut, bins, binning_field='temperature')
#print(flowline_df)
df_cut.insert(0, 'bin_center', centers[0])
#print(df_cut)
df_for_calc=df_cut[['bin_center', 'baryon_number_density', 'metallicity', 'RT_DISK_VAR_0', 'RT_DISK_VAR_1', 'RT_DISK_VAR_2', 'RT_DISK_VAR_3']]
print(df_for_calc.apply(CHF_from_array, 1))

