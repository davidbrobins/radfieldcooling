#A class to handle the usage of halo catalogs.  Includes methods to read in a halo catalog, select a subcategory of halos, select data from a halo
import yt #import yt
import pandas as pd #pandas for dataframe handling

class HaloCatalog():
    '''
    From h5 file, create a yt halo catalog and dataframe of all the data in h5filename, select between arbitrary cutoffs in some variable
    Usage:
    hc = HaloCatalog(h5filename) #Note that you may need to do this multiple times as each scale factor has 4 rockstar halo catalogs
    massive_halos=hc.select_halos(property='Mvir', min=1e9, max=1.05e9) #Cut virial mass within the desired range
    '''
    def __init__(self, a_exp): #Read in the h5 halo catalogs for the given scale factor expression (a string)
        halos_df = pd.DataFrame(columns=["ID", "Mvir", "Rvir", "x", "y", "z"]) #Create an empty dataframe for the halo catalogs
        for quadrant_index in range(4): #Loop through all 4 quadrants (indices 0,1,2,3)
            hl=yt.load("scripts/halo_catalogs/RF_"+a_exp+"_catalog/quadrant"+str(quadrant_index)+"/quadrant"+str(quadrant_index)+".0.h5") 
            halos=hl.all_data() #Extract the data from hl
            halos_df=halos_df.append(pd.DataFrame({"ID":halos["halos", "particle_identifier"], "Mvir":halos["halos", "particle_mass"].in_units("Msun/h"),'Rvir':halos["halos", "virial_radius"].in_units("kpccm/h"), "x":halos["halos", "particle_position_x"].in_units("Mpccm/h"), "y":halos['halos', "particle_position_y"].in_units("Mpccm/h"), "z":halos['halos', "particle_position_z"].in_units("Mpccm/h")}), ignore_index=True) #put halo ID,virial mass, virial radius, position in a dataframe
        self.halos_df = halos_df.sort_values("Mvir", ascending=False) #Sort values by virial mass in descending order
        self.halos_df.drop_duplicates(inplace=True)
    def select_halos(self, property='Mvir', min=1e9, max=1.1e9): #Select halos where property is between min and max
        selected_halos=self.halos_df.loc[(self.halos_df[property] >= min) & (self.halos_df[property] <= max)]
        return selected_halos

#Testing
'''
aexp='0.1451'
hc=HaloCatalog(aexp)
print(hc.halos_df)        
'''
