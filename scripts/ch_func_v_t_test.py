#Plot cooling/heating functions vs. binned T values from ch_v_t
import yt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from yt.analysis_modules.halo_analysis.api import *
import CF
h=0.6814000010490417
a=0.1451
ds = yt.load("/data/dbrobins/20/A/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art")
hl=yt.load("/home/dbrobins/repos/radfieldcooling/scripts/halo_catalogs/RF_0.1451_catalog/RF_0.1451_catalog.0.h5") #read in corresponding halo catalog
halos=hl.all_data() #extract halo information                                                                                                                                                     
halos_df=pd.DataFrame({"ID":halos["halos", "particle_identifier"], "Mvir":halos["halos", "particle_mass"].in_units("Msun/h"),'Rvir':halos["halos", "virial_radius"].in_units("kpccm/h"), "x":halos["halos", "particle_position_x"].in_units("Mpccm/h"), "y":halos['halos', "particle_position_y"].in_units("Mpccm/h"), "z":halos['halos', "particle_position_z"].in_units("Mpccm/h")}) #put halo ID,virial mass, virial radius, position in a dataframe                                                                                                                                              
halos_df=halos_df.sort_values("Mvir", ascending=False) #sort halos by virial mass in descending order
#ds.add_field(('gas', 'baryon_number_density'), function=rho_to_n_b, units='1/cm**3') #Create derived fields using functions from derived_field_ch_nb                                             
halos_mass_10e9_df=halos_df.loc[(halos_df['Mvir'] >= 10**9) & (halos_df['Mvir'] <= 1.1*10**9)] 
#ds.add_field(('gas', 'cooling_rate'), function=cooling_rate, units='erg*cm**3/s')
#ds.add_field(('gas', 'heating_rate'), function=heating_rate, units='erg*cm**3/s')
from derived_fields_ch_nb import *
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
    df = halo_sphere.to_dataframe(["baryon_number_density", "temperature","metallicity", ('artio', 'RT_DISK_VAR_0'), ('artio', 'RT_DISK_VAR_1'), ('artio', 'RT_DISK_VAR_2'), ('artio', 'RT_DISK_VAR_3'), 'cooling_rate', 'heating_rate']) #Convert halo n_b, T, Z, photoionization rates data to a pandas dataframe                 
    df = df.loc[(df['baryon_number_density']>min_den)&(df['baryon_number_density']<max_den)]  #Extract just values in the desired density bin                                                     
    df = df.sort_values(by="temperature") #Sort values by temperature (ascending)                                                                                                                 
    T=df.get('temperature').to_numpy() #Extract the three columns as numpy arrays for use later                                                                                                   
    n_b=df.get('baryon_number_density').to_numpy()
    Z=df.get('metallicity').to_numpy()
    P_LW=df.get('RT_DISK_VAR_0').to_numpy()
    P_HI=df.get('RT_DISK_VAR_1').to_numpy()
    P_HeI=df.get('RT_DISK_VAR_2').to_numpy()
    P_CVI=df.get('RT_DISK_VAR_3').to_numpy()
    cooling=df.get('cooling_rate').to_numpy()
    heating=df.get('heating_rate').to_numpy()
    num_bins=np.int(np.log10(T[-1]/T[0]))*bins_per_dex #Find number of bins by multiplying bins_per_dex by number of decades rounded to nearest integer                                           
    t_bins=np.logspace(np.log10(T[0]), np.log10(T[-1]), num_bins) #Create log-spaced bins across all temp values                                                                                  
    all_centers=10**((np.log10(t_bins[1:])+np.log10(t_bins[:-1]))/2) #Find the centers of those bins for plotting purposes
    centers=[]
    binned_c_func=[] #Arrays for cooling/heating medians in each bin                                                                                                                              
    binned_c_val=[]
    binned_h_func=[]
    binned_h_val=[]
    c_25_func=[] #Arrays for given percentile of cooling/heating in each bin                                                                                                                    
    c_25_val=[]
    h_25_func=[]
    h_25_val=[]
    c_75_func=[]
    c_75_val=[]
    h_75_func=[]
    h_75_val=[]
    cooling_stats_func=[] #Arrays to use to store all values in a given bin to then extract median/percentiles                                                                                    
    cooling_stats_val=[]
    heating_stats_func=[]
    heating_stats_val=[]
    counter = 0
    for i in range(len(t_bins)-1): #Loop through all temperature bins (more specifically, the lower limit of each bin)                                                                            
        for j in range(len(n_b)): #Loop through all n_b values (which loops through all cells)
            (cfun,hfun,ierr)=CF.frtgetcf(all_centers[i], n_b[j], Z[j], P_LW[j], P_HI[j], P_HeI[j], P_CVI[j]) #call f2py function to get gamma, lambda for the T at the center of the bin, and for the values of the other variable in the given cell
            cooling_stats_func.append(cfun)
            heating_stats_func.append(hfun)
            if t_bins[i]<=T[j] and T[j]<=t_bins[i+1]: #Check if given T value is in the bin used for the loop                                                                                     
                cooling_stats_val.append(cooling[j]) #If so, append cooling/heating values to stats arrays and add 1 to the counter                                                               
                heating_stats_val.append(heating[j])
                counter += 1
        binned_c_func.append(np.median(cooling_stats_func)) #Append the median and percentile values to the right array
        binned_h_func.append(np.median(heating_stats_func))
        c_25_func.append(np.percentile(cooling_stats_func, 25))
        h_25_func.append(np.percentile(heating_stats_func, 25))
        c_75_func.append(np.percentile(cooling_stats_func, 75))
        h_75_func.append(np.percentile(heating_stats_func, 75))
        cooling_stats_func=[] #Clear the arrays for statistics
        heating_stats_func=[]
        if counter > 0: #If the bin is non-empty, add center to used centers array, the stats values to the revlevant arrays                                                                      
            binned_c_val.append(np.median(cooling_stats_val))
            binned_h_val.append(np.median(heating_stats_val))
            c_25_val.append(np.percentile(cooling_stats_val, 25))
            h_25_val.append(np.percentile(heating_stats_val, 25))
            c_75_val.append(np.percentile(cooling_stats_val, 75))
            h_75_val.append(np.percentile(heating_stats_val, 75))
            centers.append(all_centers[i])
        counter = 0
        cooling_stats_val=[]
        heating_stats_val=[]
    ax_c.plot(all_centers, binned_c_func, label="Mvir=%.2E [M$_\odot$/h], function" % Mvir, color=colors[k]) #Plot median cooling vs. bin centers using color allocated to that halo, label with mass to 3 sf 
    #ax_c.fill_between(all_centers, c_25_func, c_75_func, alpha=0.2, color=colors[k]) #Fill between 25th-75th percentiles in same color                                                           
    ax_h.plot(all_centers, binned_h_func, color=colors[k]) #Do the same thing for heating rate                                                                                                    
    #ax_h.fill_between(all_centers, h_25_func, h_75_func, alpha=0.2, color=colors[k])
    ax_c.plot(centers, binned_c_val, label="Mvir=%.2E [M$_\odot$/h], values" % Mvir, color=colors[k], linestyle='dashed') #Plot median cooling vs. bin centers using color allocated to that halo, label with mass to 3 sig. figs 
    #ax_c.fill_between(centers, c_25, c_75, alpha=0.2, color=colors[k]) #Fill between 25th-75th percentiles in same color                                                                        
    ax_h.plot(centers, binned_h_val, color=colors[k], linestyle='dashed') #Do the same thing for heating rate                                                                                    
    #ax_h.fill_between(centers, h_25_val, h_75_val, alpha=0.2, color=colors[k])
    ax_c.set_xscale('log') #Set all axes to log scale                                                                                                                                             
    ax_c.set_yscale('log')
    ax_h.set_xscale('log')
    ax_h.set_yscale('log')
ax_c.set_xlabel("Temperature T [K]") #Label the axes, create a legend, save the plot                                                                                                              
ax_h.set_xlabel("Temperature T [K]")
ax_c.set_ylabel("$\Lambda$ [erg cm$^3$/s]")
ax_h.set_ylabel("$\Gamma$ [erg cm$^3$/s]")
fig.legend(prop={'size':5})
fig.savefig("/home/dbrobins/repos/radfieldcooling/plots/ch_comp_v_t_mass_10e9.png")
