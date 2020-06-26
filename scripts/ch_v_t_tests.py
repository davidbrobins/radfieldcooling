#Plots temperature-binned values of cooling/heating rates evaluated for each cell
import yt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from yt.analysis_modules.halo_analysis.api import *
from derived_fields_ch_nb import *
h=0.6814000010490417
a=0.1451
ds = yt.load("/data/dbrobins/20/A/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art")
hl=yt.load("/home/dbrobins/repos/radfieldcooling/scripts/halo_catalogs/RF_0.1451_catalog/RF_0.1451_catalog.0.h5") #read in corresponding halo catalog 
halos=hl.all_data() #extract halo information
halos_df=pd.DataFrame({"ID":halos["halos", "particle_identifier"], "Mvir":halos["halos", "particle_mass"].in_units("Msun/h"),'Rvir':halos["halos", "virial_radius"].in_units("kpccm/h"), "x":halos["halos", "particle_position_x"].in_units("Mpccm/h"), "y":halos['halos', "particle_position_y"].in_units("Mpccm/h"), "z":halos['halos', "particle_position_z"].in_units("Mpccm/h")}) #put halo ID, virial mass, virial radius, position in a dataframe
halos_df=halos_df.sort_values("Mvir", ascending=False) #sort halos by virial mass in descending order
#halos_mass_10e9_df=halos_df.loc[(halos_df['Mvir'] >= 10**9) & (halos_df['Mvir'] <= 1.1*10**9)] 
ds.add_field(('gas', 'baryon_number_density'), function=rho_to_n_b, units='1/cm**3') #Create derived fields using functions from derived_field_ch_nb
ds.add_field(('gas', 'cooling_rate'), function=cooling_rate, units='erg*cm**3/s')
ds.add_field(('gas', 'heating_rate'), function=heating_rate, units='erg*cm**3/s')
min_den=1
max_den=10 #min and max density in 1/cm^3
bins_per_dex=20 #number of bins per factor of 10 (dex) in T (K)
fig, (ax_c, ax_h)=plt.subplots(2,1,sharex=True) #Create two horizontal plots with the same x-axis (temperature) for the cooling and heating rates
colors = ['red', 'blue', 'orange'] #Create an array of colors
for k in range(3):
    ID=halos_df.iloc[k]['ID'] #extract halo ID, virial mass, virial radius, halo position                                                                                               
    Mvir=halos_df.iloc[k]['Mvir']
    R=halos_df.iloc[k]['Rvir']
    x=halos_df.iloc[k]['x']
    y=halos_df.iloc[k]['y']
    z=halos_df.iloc[k]['z']
    halo_sphere = ds.sphere(yt.YTArray([x*a/h, y*a/h, z*a/h], "Mpc"),(R*a/h, "kpc")) #pick out test halo                                                                                      
    df = halo_sphere.to_dataframe(["baryon_number_density", "temperature","cooling_rate", "heating_rate"]) #Convert halo n_b, cooling+heating function data to a pandas dataframe             
    df = df.loc[(df['baryon_number_density']>min_den)&(df['baryon_number_density']<max_den)]  #Extract just values in the desired density bin                                                 
    df = df.sort_values(by="temperature") #Sort values by temperature (ascending)
    T=df.get('temperature').to_numpy() #Extract the three columns as numpy arrays for use later                                                                                   
    cooling=df.get('cooling_rate').to_numpy()
    heating=df.get('heating_rate').to_numpy()
    num_bins=np.int(np.log10(T[-1]/T[0]))*bins_per_dex #Find number of bins by multiplying bins_per_dex by number of decades rounded to nearest integer
    t_bins=np.logspace(np.log10(T[0]), np.log10(T[-1]), num_bins) #Create log-spaced bins across all temp values
    all_centers=10**((np.log10(t_bins[1:])+np.log10(t_bins[:-1]))/2) #Find the centers of those bins for plotting purposes
    centers = [] #Create array for centers of non-empty bins
    binned_c=[] #Arrays for cooling/heating medians in each bin
    binned_h=[]
    c_25=[] #Arrays for given percentile of cooling/heating in each bin
    h_25=[]
    c_75=[]
    h_75=[]
    cooling_stats=[] #Arrays to use to store all values in a given bin to then extract median/percentiles
    heating_stats=[]
    counter = 0 #Counter to check whether a bin is empty
    for i in range(len(t_bins)-1): #Loop through all temperature bins (more specifically, the lower limit of each bin)
        for j in range(len(T)): #Loop through all T values
            if t_bins[i]<=T[j] and T[j]<=t_bins[i+1]: #Check if given T value is in the bin used for the loop
                cooling_stats.append(cooling[j]) #If so, append cooling/heating values to stats arrays and add 1 to the counter
                heating_stats.append(heating[j])
                counter += 1
        print(all_centers[i], counter) #Print as a check
        if counter > 0: #If the bin is non-empty, add center to used centers array, the stats values to the revlevant arrays
            binned_c.append(np.median(cooling_stats))
            binned_h.append(np.median(heating_stats))
            c_25.append(np.percentile(cooling_stats, 25))
            h_25.append(np.percentile(heating_stats, 25))
            c_75.append(np.percentile(cooling_stats, 75))
            h_75.append(np.percentile(heating_stats, 75))
            centers.append(all_centers[i])
        counter = 0 #Reset the counter and the stats arrays
        cooling_stats=[]
        heating_stats=[]
    ax_c.plot(centers, binned_c, label="Mvir=%.2E [M$_\odot$/h]" % Mvir, color=colors[k]) #Plot median cooling vs. bin centers using color allocated to that halo, label with mass to 3 sig. figs 
    ax_c.fill_between(centers, c_25, c_75, alpha=0.2, color=colors[k]) #Fill between 25th-75th percentiles in same color
    ax_h.plot(centers, binned_h, color=colors[k]) #Do the same thing for heating rate
    ax_h.fill_between(centers, h_25, h_75, alpha=0.2, color=colors[k])
    ax_c.set_xscale('log') #Set all axes to log scale
    ax_c.set_yscale('log')
    ax_h.set_xscale('log')
    ax_h.set_yscale('log')
    # Now do the same for 3 halos at ~10^9 Msun/h
    '''
    ID=halos_mass_10e9_df.iloc[k]['ID'] #extract halo ID, virial mass, virial radius, halo position                                                                                                          
    Mvir=halos_mass_10e9_df.iloc[k]['Mvir']
    R=halos_mass_10e9_df.iloc[k]['Rvir']
    x=halos_mass_10e9_df.iloc[k]['x']
    y=halos_mass_10e9_df.iloc[k]['y']
    z=halos_mass_10e9_df.iloc[k]['z']
    halo_sphere = ds.sphere(yt.YTArray([x*a/h, y*a/h, z*a/h], "Mpc"),(R*a/h, "kpc")) #pick out test halo                                                                                           
    df = halo_sphere.to_dataframe(["baryon_number_density", "temperature","cooling_rate", "heating_rate"]) #Convert halo n_b, cooling+heating function data to a pandas dataframe                  
    df = df.loc[(df['baryon_number_density']>min_den)&(df['baryon_number_density']<max_den)]  #Extract just values in the desired density bin                                                      
    df = df.sort_values(by="temperature")
    T=df.get('temperature').to_numpy() #Extract the three columns as numpy arrays for use later                                                                                                    
    cooling=df.get('cooling_rate').to_numpy()
    heating=df.get('heating_rate').to_numpy()
    num_bins=np.int(np.log10(T[-1]/T[0]))*bins_per_dex
    t_bins=np.logspace(np.log10(T[0]), np.log10(T[-1]), num_bins)
    all_centers=10**((np.log10(t_bins[1:])+np.log10(t_bins[:-1]))/2)
    centers = []
    binned_c=[]
    binned_h=[]
    c_25=[]
    h_25=[]
    c_75=[]
    h_75=[]
    cooling_stats=[]
    heating_stats=[]
    counter = 0
    for i in range(len(t_bins)-1):
        for j in range(len(T)):
            if t_bins[i]<=T[j] and T[j]<=t_bins[i+1]:
                cooling_stats.append(cooling[j])
                heating_stats.append(heating[j])
                counter += 1
        print(all_centers[i], counter)
        if counter > 0:
            binned_c.append(np.median(cooling_stats))
            binned_h.append(np.median(heating_stats))
            c_25.append(np.percentile(cooling_stats, 25))
            h_25.append(np.percentile(heating_stats, 25))
            c_75.append(np.percentile(cooling_stats, 75))
            h_75.append(np.percentile(heating_stats, 75))
            centers.append(all_centers[i])
        counter = 0
        cooling_stats=[]
        heating_stats=[]
    ax_c.plot(centers, binned_c, label="Mvir=%.2E [M$_\odot$/h]" % Mvir)
    if k==1:
        ax_c.plot(centers, c_25, linestyle='dashed')
        ax_c.plot(centers, c_75, linestyle='dashed')
    ax_h.plot(centers, binned_h)
    if k==1:
        ax_h.plot(centers, h_25, linestyle='dashed')
        ax_h.plot(centers, h_75, linestyle='dashed')
    ax_c.set_xscale('log')
    ax_c.set_yscale('log')
    ax_h.set_xscale('log')
    ax_h.set_yscale('log')
    '''
ax_c.set_xlabel("Temperature T [K]") #Label the axes, create a legend, save the plot
ax_h.set_xlabel("Temperature T [K]")
ax_c.set_ylabel("$\Lambda$ [erg cm$^3$/s]")
ax_h.set_ylabel("$\Gamma$ [erg cm$^3$/s]")
fig.legend(prop={'size':5})
fig.savefig("/home/dbrobins/repos/radfieldcooling/plots/ch_v_t_most_massive.png")

