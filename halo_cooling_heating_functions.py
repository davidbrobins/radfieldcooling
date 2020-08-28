# A module to do binning and create isochronous and flowline averaged cooling/heating function arrays for a given halo
import yt #We need yt to do this
import numpy as np #Numpy
import pandas as pd #Pandas (dataframe handling)
import CF #f2py module that contains python wrappers for calculating cooling and heating functions
from halo_catalog import HaloCatalog #Module defining HaloCatalog() class
h=0.6814000010490417 #The hubble constant: H0=100*h km/s/Mpc
k_boltz=1.3807e-16 #Boltzmann's constant in erg/K
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

def get_median_percentiles(binning_field, bin_center, chft_df):
    '''
    Function to take a dataframe of cooling_rate, heating_rate, chft_df and take the median, 25th and 75th percentiles of each column
    Inputs:
    binning_field (str): The name of the binning field
    bin_center (array): The center of the relevant bin
    chft_df (dataframe): A dataframe with three columns: cooling rate, heating rate, cooling time (in that order)
    Outputs:
    output_df: A one-row dataframe with columns for binning_field (center of each bin), then median, 25th and 75th percentiles for CF, HF, and cooling time
    '''
    cf_med, hf_med, ct_med = chft_df.apply(np.median, 0) #Get the median of each column and store                                                                                       
    cf_25, hf_25, ct_25 = chft_df.apply(np.percentile, 0, args=(25,)) #Get the 25th percentile of each column                                                                            
    cf_75, hf_75, ct_75 = chft_df.apply(np.percentile, 0, args=(75,)) #Get the 75th percentile of each column                                                                            
    averages_df={binning_field:bin_center, 'CF_median':cf_med, 'CF_25':cf_25, 'CF_75':cf_75, 'HF_median':hf_med, 'HF_25':hf_25, 'HF_75':hf_75, 'ctime_m\
edian':ct_med, 'ctime_25':ct_25, 'ctime_75':ct_75} #Put the resulting values in an output dataframe
    return averages_df
       
def flowline_averaging(df, bins, binning_field='temperature', to_save=True, filepath='saved_dataframes/a0.1451/newest'): 
    '''
    Function to take a dataframe and given bins for binning_field and return flowline averages (and 25th, 75th percentiles) for cooling_rate, heating_rate, cooling_time in those bins
    Inputs:
    df (dataframe): The dataframe used for averaging, must (at minimum) contain columns of binning_field, cooling_rate, heating_rate, and cooling_time
    bins (ndarray): Array of (logspaced) bins including left and right edges
    binning_field (string): The field to bin in (default to temperature)
    to_save (bool): Whether or not to save the resulting arrays (defaults to True)
    filepath (str): A string giving a filepath and name for the saved pickle file
    Outputs: 
    flowline_averages_df (dataframe): dataframe with centers of bin in binning_field and median, 25th, 75th percentile for cooling rate, heating rate, and cooling time
    '''
    centers=get_bin_centers(bins) #Get (log-spaced) centers of the bins
    flowline_averages_df=pd.DataFrame(columns=[binning_field, 'CF_median', 'CF_25', 'CF_75', 'HF_median', 'HF_25', 'HF_75', 'ctime_median', 'ctime_25', 'ctime_75'])
    for bin_index in range(len(bins)-1): #Loop through the left edges of each bin
        binned_df=cut_fields_data(df, property=binning_field, min=bins[bin_index], max=bins[bin_index+1]) #Cut df to only contain the data where binning_field is within the appropriate bin
        cols_for_averaging=binned_df[['cooling_rate', 'heating_rate', 'cooling_time']]
        if len(cols_for_averaging.index)>0: #If cols_for_averaging has at least 1 row
            stats_df=get_median_percentiles(binning_field, centers[bin_index], cols_for_averaging) #Extract stats
            flowline_averages_df=flowline_averages_df.append(stats_df, ignore_index=True)
        else:
            flowline_averages_df=flowline_averages_df.append({binning_field:centers[bin_index], 'CF_median':np.nan, 'CF_25':np.nan, 'CF_75':np.nan, 'HF_median':np.nan, 'HF_25':np.nan, 'HF_75':np.nan, 'ctime_median':np.nan, 'ctime_25':np.nan, 'ctime_75':np.nan}, ignore_index=True)
    if to_save==True:
        flowline_averages_df.to_pickle(filepath+'_flowline.pkl')
    return flowline_averages_df #Return the dataframe after the loop completes
    #if savenpz==True: save centers, c_func_median, h_func_median, c_time_median, percentiles

def get_CHFT(frtargs, cols_to_use=['T', 'n_b', 'Z', 'P_LW', 'P_HI', 'P_HeI', 'P_CVI']): #get cooling function from dataframe of input values (for use in calculations for isochronous binning)          
    '''
    Function to get cooling and heating functions from an input array
    Input:
    frtargs (dataframe/dataframe row): Row containing [temperature, baryon number density, metallicity, P_LW (RT_DISK_VAR_0), P_HI (RT_DISK_VAR_1), P_HeI (RT_DISK_VAR_2), P_CVI (RT_DISK_VAR_3)]
    cols_to_use (array): Names of the columns to use, can vary depending on setup (but default to 'T', 'n_b', 'Z', 'P_LW', 'P_HI', 'P_HeI', 'P_CVI')
    Outputs:
    cfun, hfun, c_time (floats): The cooling and heating functions, and the cooling time at those values
    '''
    cfargs=[frtargs[col] for col in cols_to_use]
    (cfun,hfun,ierr)=CF.frtgetcf(*cfargs) #call f2py function to get cooling and heating functions for these inputs
    c_time = k_boltz*cfargs[0]/(cfargs[1]*cfun) #(kB*T/n_b*Lambda(T))
    return [cfun, hfun, c_time] #Return the cooling and heating functions, and the cooling time

def isochronous_averaging(df, bins, binning_field='temperature', to_save=True, filepath='saved_dataframes/a0.1451/newest'): #add a boolean argument to_save that defaults to True
    '''
    Function to take a dataframe and given bins for binning_field and return isochronous averages (and 25th, 75th percentiles) for cooling_rate, heating_rate, cooling_time in those bins                  
    Inputs:                                                                                                                                                                              
    df (dataframe): The dataframe used for averaging, must (at minimum) contain columns of T, n_b, Z, P_LW, P_HI, P_HeI, P_CVI                   
    bins (ndarray): Array of (logspaced) bins including left and right edges                                                                                                                               
    binning_field (string): The field to bin in (default to temperature), should be either temperature or baryon_number_density                                                                      
    to_save (bool): Whether or not to save the resulting arrays (defaults to True)                                                                                                                                                                                                                  
    filepath (str): A string giving a filepath and name for the saved pickle file                                                                                                                     
    Outputs:                                                                                                                                                                                               
    isochronous_averages_df (dataframe): dataframe with centers of bin in binning_field and median, 25th, 75th percentile for cooling rate, heating rate, and cooling time                                 
    ''' 
    centers=get_bin_centers(bins) #Get bin centers
    isochronous_averages_df=pd.DataFrame(columns=[binning_field, 'CF_median', 'CF_25', 'CF_75', 'HF_median', 'HF_25', 'HF_75', 'ctime_median', 'ctime_25', 'ctime_75'])   
    for bin_center in centers: #Loop through the center of each bin
        df['bin_center']=bin_center #Add a column to the dataframe with bin_center as the constant value of the column
        if binning_field=='temperature': #Array of keys must be [T, n_b, Z, P_LW, P_HI, P_HeI, P_CVI].  Set binning_field to 'bin_center', others to appropriate yt field
            args_for_chft=['bin_center', 'baryon_number_density', 'metallicity', ('artio','RT_DISK_VAR_0'), ('artio','RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio','RT_DISK_VAR_3')]
        elif binning_field=='baryon_number_density':
            args_for_chft=['temperature', 'bin_center', 'metallicity', ('artio','RT_DISK_VAR_0'), ('artio,','RT_DISK_VAR_1'), ('artio','RT_DISK_VAR_2'), ('artio','RT_DISK_VAR_3')]
        df_for_stats=df.apply(get_CHFT, 1, cols_to_use=args_for_chft, result_type='expand') #Apply get_CHFT to each row and expand to a dataframe
        binned_stats_df=get_median_percentiles(binning_field, bin_center, df_for_stats)
        isochronous_averages_df=isochronous_averages_df.append(binned_stats_df, ignore_index=True)
    if to_save==True:
        isochronous_averages_df.to_pickle(filepath+'_isochronous.pkl')
    return isochronous_averages_df
    #binning_field should be either T ('temperature') or n_b  ('baryon_number_density')
    #loop through centers, add a column to df which is just that T value
    #to df['T_for_calc', 'baryon_number_density', 'metallicity', 'RT_DISK_VAR_0', 'RT_DISK_VAR_1', 'RT_DISK_VAR_2', 'RT_DISK_VAR_3'].apply()
    #Take median of the result
    #if savenpz==True: save the results
    

#testing:
'''
import matplotlib.pyplot as plt
fig, (ax_c, ax_h)=plt.subplots(2,1,sharex=True) #Create two horizontal plots with the same x-axis (temperature) for the cooling and heating rates                                                          
colors = ['red', 'blue', 'orange'] #Create an array of colors 
hc=HaloCatalog("/home/dbrobins/repos/radfieldcooling/radfieldcooling/scripts/halo_catalogs/RF_0.1451_catalog/RF_0.1451_catalog.0.h5")
halos_df=hc.halos_df
for massive_halo in range(3):
    Mvir=halos_df.iloc[massive_halo]['Mvir']
    df = get_halo_sphere_data(halos_df, "/nfs/turbo/lsa-cavestru/shared_data/CROC/dbrobins/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art", massive_halo, 0.1451, ['temperature', 'baryon_number_density','metallicity', ('artio', 'RT_DISK_VAR_0'), ('artio', 'RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio', 'RT_DISK_VAR_3'), 'cooling_rate', 'heating_rate', 'cooling_time'])
    df_cut = cut_fields_data(df, property='baryon_number_density', min=1, max=10) 
    df_cut = df_cut.sort_values(by='temperature')
    bins=log_spaced_bins(df_cut['temperature'].to_numpy()[0], df_cut['temperature'].to_numpy()[-1], 20)
    flowline_df = flowline_averaging(df_cut, bins, binning_field='temperature', to_save=True, name='most_massive_'+str(massive_halo))
    print(flowline_df)
    isochronous_df=isochronous_averaging(df_cut, bins, to_save=True, name='most_massive_'+str(massive_halo))
    print(isochronous_df)
    ax_c.plot(isochronous_df['temperature'], isochronous_df['CF_median'], label="Mvir=%.2E [M$_\odot$/h], isochronous" % Mvir, color=colors[massive_halo]) #Plot median cooling vs. bin centers using color allocated to that halo, label with mass to 3 sf                                                                                                                                                           
    ax_h.plot(isochronous_df['temperature'], isochronous_df['HF_median'], color=colors[massive_halo]) #Do the same thing for heating rate                                                                  
    ax_c.plot(flowline_df['temperature'], flowline_df['CF_median'], label="Mvir=%.2E [M$_\odot$/h], values" % Mvir, color=colors[massive_halo], linestyle='dashed') #Plot median cooling vs. bin centers using color allocated to that halo, label with mass to 3 sig. figs                                                                                                                                           
    ax_h.plot(flowline_df['temperature'], flowline_df['HF_median'], color=colors[massive_halo], linestyle='dashed') #Do the same thing for heating rate                                                    
    ax_c.set_xscale('log') #Set all axes to log scale                                                                                                                                                      
    ax_c.set_yscale('log')
    ax_h.set_xscale('log')
    ax_h.set_yscale('log')
ax_c.set_xlabel("Temperature T [K]") #Label the axes, create a legend, save the plot                                                                                                                      
ax_h.set_xlabel("Temperature T [K]")
ax_c.set_ylabel("$\Lambda$ [erg cm$^3$/s]")
ax_h.set_ylabel("$\Gamma$ [erg cm$^3$/s]")
fig.legend(prop={'size':5})
fig.savefig("/home/dbrobins/repos/radfieldcooling/plots/refactored_flowline_isochronous_comparison_T.png")
'''
