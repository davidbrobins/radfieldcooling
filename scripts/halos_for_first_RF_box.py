# Code to create a yt h5  halo catalog for the radiation field box at a=0.1451
# Only needs to be run once (it's already been run)
import yt
from yt.analysis_modules.halo_analysis.api import *
ds = yt.load("/data/dbrobins/20/A/rei20c1_a0.1451_RF/rei20c1_a0.1451_RF.art") #Read in the whole box
halos_ds=yt.load("/data/dbrobins/20/A/rs/halos_10.0.bin") #Read in the binary halo catalog at the right scale factor
hc=HaloCatalog(data_ds=ds, halos_ds=halos_ds, output_dir="halo_catalogs/RF_0.1451_catalog") #Make the halo catalog and set the directory for the h5 file
hc.create() #Create the h5 file
