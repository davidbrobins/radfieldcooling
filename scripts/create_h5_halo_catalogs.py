# Code to create  yt h5  halo catalog for all the radiation field boxes (4 for each box corresponding to the random quadrants)
# Only needs to be run once (it's already been run)
import yt
from yt.analysis_modules.halo_analysis.api import *
a_dict={'0.0905':'3', '0.1112':'5', '0.1284':'7', '0.1451':'10', '0.1667':'14'} #Dictionary with entries of the form a_exp:index_exp
for a_exp in a_dict.keys():
    ds=yt.load("/nfs/turbo/lsa-cavestru/shared_data/CROC/dbrobins/rei20c1_a"+a_exp+"_RF/rei20c1_a"+a_exp+"_RF.art") #Read in the whole box
    for quadrant_index in range(4):
        halos_ds=yt.load("/nfs/turbo/lsa-cavestru/shared_data/CROC/dbrobins/halo_catalogs/halos_"+a_dict[a_exp]+"."+str(quadrant_index)+".bin") #Read in the halo catalogs for each quadrant
        hc=HaloCatalog(data_ds=ds, halos_ds=halos_ds, output_dir="halo_catalogs/RF_"+a_exp+"_catalog/quadrant"+str(quadrant_index)) #Send the halo catalog to the right directory
        hc.create() #Create the h5 file
