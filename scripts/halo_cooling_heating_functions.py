# A module to do binning and create isochronous and flowline averaged cooling/heating function arrays for a given halo
import yt #We need yt to do this
import numpy as np #Numpy
import pandas as pd #Pandas (dataframe handling)
import CF #f2py module that contains python wrappers for calculating cooling and heating functions
import halo_catalog #Module defining HaloCatalog() class
from derived_fields_ch_nb import * #Import module with baryon_number_density, cooling_rate, heating_rate, cooling_time derived fields

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

def get_halo_sphere_data(halos_df, ds, halo_index, a, h, fields):
    '''
    Function to extract a dataframe with given fields and halo info dictionary for a given halo index
    Inputs:
    halos_df (dataframe): A dataframe with the given halos, likely imported from HaloCatalog() class
    ds (ytObject): The yt data set being used
    halo_index (int): The index of the halo in halos_df 
    a (float): scale factor [might be removed]
    h (float): Hubble parameter (dimensionless) [might be removed]
    fields (array of strings): Fields to be selected.  Choose from 'baryon_number_density', 'temperature', 'metallicity', ('artio', 'RT_DISK_VAR_[0-3]'), 'cooling_rate', 'heating_rate', 'cooling_time'
    Outputs:
    df (dataframe): Dataframe containing fields as the column labels
    '''
    halo_info=halos_df.iloc[halo_index] #Extract a dictionary with halo ID, Mvir, Rvir, x, y, z
    halo_sphere=ds.sphere(yt.YTArray([halo_info['x']*a/h, halo_info['y']*a/h, halo_info['z']*a/h], "Mpc", (halo_info['Rvir']*a/h, "kpc")) #Pick out a Rvir-radius sphere centered on the halo center
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

