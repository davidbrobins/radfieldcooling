#A class to handle the usage of halo catalogs.  Includes methods to read in a halo catalog, select a subcategory of halos, select data from a halo
#Data: T, n_b, Z, P_LW, P_HI, P_HeI, P_CVI, cooling function, heating function, cooling time (can add more fields if it becomes necessary  
import yt #import yt
import pandas as pd #pandas for dataframe handling

class HaloCatalog():
    '''
    From h5 file, create a yt halo catalog and dataframe of all the data in h5filename, select between arbitrary cutoffs in some variable
    Usage:
    hc = HaloCatalog(h5filename) #Note that you may need to do this multiple times as each scale factor has 4 rockstar halo catalogs
    massive_halos=hc.select_halos(property='Mvir', min=1e9, max=1.05e9) #Cut virial mass within the desired range
    '''
    def __init__(self, h5filename): #Read in the h5 halo catalog h5filename (a string) and convert to a dataframe
        hl=yt.load(h5filename) #h5filename should be in form "/home/dbrobins/repos/radfieldcooling/scripts/halo_catalogs/RF_"+aexp+"_catalog/RF_"+aexp+"_catalog."+index+".h5", where index is 0-3
        halos=hl.all_data() #Extract the data from hl
        halos_df = pd.DataFrame({"ID":halos["halos", "particle_identifier"], "Mvir":halos["halos", "particle_mass"].in_units("Msun/h"),'Rvir':halos["halos", "virial_radius"].in_units("kpccm/h"), "x":halos["halos", "particle_position_x"].in_units("Mpccm/h"), "y":halos['halos', "particle_position_y"].in_units("Mpccm/h"), "z":halos['halos', "particle_position_z"].in_units("Mpccm/h")}) #put halo ID,virial mass, virial radius, position in a dataframe
        self.halos_df = halos_df.sort_values("Mvir", ascending=False) #Sort values by virial mass in descending order

    def select_halos(self, property='Mvir', min=1e9, max=1.1e9): #Select halos where property is between min and max
        selected_halos=self.halos_df.loc[(self.halos_df[property] >= min) & (self.halos_df[property] <= max)]
        return selected_halos

#Testing
#aexp='0.1451'
#hc=HaloCatalog("/home/dbrobins/repos/radfieldcooling/scripts/halo_catalogs/RF_"+aexp+"_catalog/RF_"+aexp+"_catalog.0.h5")
#print(hc.select_halos(property='Rvir', min=1, max=25))        
