#Plots cooling/heating median/spread in 8 density bins per dex in 3 halos in most massive, 10^9 Msun/h bins
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
halos_mass_10e9_df=halos_df.loc[(halos_df['Mvir'] >= 10**9) & (halos_df['Mvir'] <= 1.1*10**9)]
ds.add_field(('gas', 'baryon_number_density'), function=rho_to_n_b, units='1/cm**3') #Add derived fields using functions from derived_fields_ch_nb
ds.add_field(('gas', 'cooling_rate'), function=cooling_rate, units='erg*cm**3/s')
ds.add_field(('gas', 'heating_rate'), function=heating_rate, units='erg*cm**3/s')
for k in range(3):
    ID=halos_df.iloc[k]['ID'] #extract halo ID, virial mass, virial radius, halo position
    Mvir=halos_df.iloc[k]['Mvir']
    R=halos_df.iloc[k]['Rvir']
    x=halos_df.iloc[k]['x']
    y=halos_df.iloc[k]['y']
    z=halos_df.iloc[k]['z']
    halo_sphere = ds.sphere(yt.YTArray([x*a/h, y*a/h, z*a/h], "Mpc"),(R*a/h, "kpc")) #pick out test halo  
    df = halo_sphere.to_dataframe(["baryon_number_density", "cooling_rate", "heating_rate"]) #Convert halo n_b, cooling+heating function data to a pandas dataframe
    df = df.sort_values(by='baryon_number_density') #Sort the dataframe by n_b
    n_b=df.get('baryon_number_density').to_numpy() #Extract the three columns as numpy arrays for easier use later
    cooling=df.get('cooling_rate').to_numpy() 
    heating=df.get('heating_rate').to_numpy()
    num_bins= np.int(np.log10(n_b[-1]/n_b[0]))*8 #Set the number of density bins as 15 times the number of powers of 10 spanned by the min and max n_b
    bins=np.logspace(np.log10(n_b[0]), np.log10(n_b[-1]), num=num_bins) #Create that number of log-spaced density bins between min and max n_b
    centers = 10**( (np.log10(bins[1:]) + np.log10(bins[:-1])) / 2 ) #Find the center of each density bin
    binned_cooling=[] #Create array to hold the values of median, spread for each density bin
    binned_heating=[]
    cooling_spread=[]
    heating_spread=[]
    cooling_for_stats=[] #Create arrays to hold cooling, heating func values to get median, spread in each bin
    heating_for_stats=[]
    binned_density=[]
    for i in range(len(bins)-1): #Loop through the density bins
        for j in range(len(n_b)): #Loop through all density values
            if bins[i] <= n_b[j] and n_b[j] <= bins[i+1]: #Check if the given density value is inside the given density bin
                cooling_for_stats.append(cooling[j]) #If so, append the cooling and heating function values to the arrays for that bin
                heating_for_stats.append(heating[j])
        binned_cooling.append(np.median(cooling_for_stats)) #Compute the median, spread in cooling and heating and append to the relevant arrays
        binned_heating.append(np.median(heating_for_stats))
        cooling_spread.append(np.percentile(cooling_for_stats, 75)-np.percentile(cooling_for_stats, 25))
        heating_spread.append(np.percentile(heating_for_stats, 75)-np.percentile(heating_for_stats, 25))
        cooling_for_stats=[] #Clear the arrays used to calculate the statistics in each bin
        heating_for_stats=[]
    directory="/home/dbrobins/repos/radfieldcooling/plots/cooling_heating_prelim_tests_5_22/" #The path to store plots in
    plt.plot(centers, heating_spread, label="Mvir: %e [M$_\odot$/h]" % Mvir, color='green') #Plot spread on log-log plot, labeling with virial mass
    plt.xscale('log')
    plt.yscale('log')
# Do it all again for three most massive halos in 10^9 Msun/h mass bin
for k in range(3):
    ID=halos_mass_10e9_df.iloc[k]['ID'] #extract halo ID, virial mass, virial radius, halo position                                                                                                  
    Mvir=halos_mass_10e9_df.iloc[k]['Mvir']
    R=halos_mass_10e9_df.iloc[k]['Rvir']
    x=halos_mass_10e9_df.iloc[k]['x']
    y=halos_mass_10e9_df.iloc[k]['y']
    z=halos_mass_10e9_df.iloc[k]['z']
    halo_sphere = ds.sphere(yt.YTArray([x*a/h, y*a/h, z*a/h], "Mpc"),(R*a/h, "kpc")) #pick out test halo                                                                                             
    df = halo_sphere.to_dataframe(["baryon_number_density", "cooling_rate", "heating_rate"]) #Convert halo n_b, cooling+heating function data to a pandas dataframe                                  
    df = df.sort_values(by='baryon_number_density') #Sort the dataframe by n_b                                                                                                                       
    n_b=df.get('baryon_number_density').to_numpy() #Extract the three columns as numpy arrays for easier use later                                                                                   
    cooling=df.get('cooling_rate').to_numpy()
    heating=df.get('heating_rate').to_numpy()
    num_bins= np.int(np.log10(n_b[-1]/n_b[0]))*8 #Set the number of density bins as 15 times the number of powers of 10 spanned by the min and max n_b                                            
    bins=np.logspace(np.log10(n_b[0]), np.log10(n_b[-1]), num=num_bins) #Create that number of log-spaced density bins between min and max n_b                                                       
    centers = 10**( (np.log10(bins[1:]) + np.log10(bins[:-1])) / 2 )
    binned_cooling=[]
    binned_heating=[]
    cooling_spread=[]
    heating_spread=[]
    cooling_for_stats=[]
    heating_for_stats=[]
    binned_density=[]
    for i in range(len(bins)-1):
        for j in range(len(n_b)):
            if bins[i] <= n_b[j] and n_b[j] <= bins[i+1]:
                cooling_for_stats.append(cooling[j])
                heating_for_stats.append(heating[j])
        #print(cooling_for_stats)
        binned_cooling.append(np.median(cooling_for_stats))
        binned_heating.append(np.median(heating_for_stats))
        cooling_spread.append(np.percentile(cooling_for_stats, 75)-np.percentile(cooling_for_stats, 25))
        heating_spread.append(np.percentile(heating_for_stats, 75)-np.percentile(heating_for_stats, 25))
        cooling_for_stats=[]
        heating_for_stats=[]
    plt.plot(centers, heating_spread, label="Mvir: %e [M$_\odot$/h]" % Mvir, color='blue')
    plt.xscale('log')
    plt.yscale('log')
plt.xlabel("Baryon number density [1/cm$^3$]") #Label the axes of the plot
plt.ylabel("Spread in heating rate [erg cm$^3$/s]") 
plt.ylim([1e-29, 1e-21]) #Set limits on y axis to make it easier to view (found by trial and error)
plt.legend() #Create a legend
plt.savefig(directory+"h_spread_comparison") #Save the plot with a relevant title
    
    

