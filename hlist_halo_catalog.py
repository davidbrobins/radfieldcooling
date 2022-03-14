#A class to handle the usage of the rockstar halo catalogs in hlists.  Includes methods to read in a halo catalog, select a subcategory of halos, select data from a halo

import pandas as pd #pandas for dataframe handling

aexp_to_hlist={'0.0905':'0.09049', '0.1002':'0.10017', '0.1112':'0.11121', '0.1284':'0.12844', '0.1400':'0.14005', '0.1451':'0.14512', '0.1505':'0.15047', '0.1552':'0.15515', '0.1601':'0.16005', '0.1667':'0.16671'}

class hlistHaloCatalog():
    '''
    From .list file, create a dataframe of all the data in the .list file, select between arbitrary cutoffs in some variable
    Usage:
    hc = HaloCatalog(a_exp) #Note that you may need to do this multiple times as each scale factor has 4 rockstar halo catalogs
    massive_halos=hc.select_halos(property='Mvir', min=1e9, max=1.05e9) #Cut virial mass within the desired range
    '''
    def __init__(self, a_exp, halo_catalog_dir='/nfs/turbo/lsa-cavestru/shared_data/CROC/dbrobins/hlist_halo_catalogs/'): #Read in the hlist halo catalog for the given scale factor expression (a string)
        #put scale factor; halo ID; parent ID (-1 if central halo, >0 for subhalo); virial mass (Msun/h); virial radius (comoving kpc/h); x,y,z position (comoving Mpc/h) in a dataframe
        halos_df = pd.read_csv(halo_catalog_dir+'hlist_'+aexp_to_hlist[a_exp]+'.list', skiprows=64, names=['a', 'ID', 'pid', 'Mvir', 'Rvir', 'x', 'y', 'z'], delim_whitespace=True, usecols=[0, 1, 5, 10, 11, 17, 18, 19], index_col=False)
        self.halos_df = halos_df.sort_values("Mvir", ascending=False) #Sort values by virial mass in descending order
    def select_halos(self, property='Mvir', min=1e9, max=9.9e9): #Select halos where property is between min and max
        selected_halos=self.halos_df.loc[(self.halos_df[property] >= min) & (self.halos_df[property] <= max)]
        return selected_halos

#Testing
'''
aexp='0.0905'
hc=hlistHaloCatalog(aexp)
halos_df=hc.halos_df
halos_df=halos_df.loc[halos_df['pid']>-1]
print(halos_df)
'''
