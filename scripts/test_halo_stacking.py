#Script to test whether stacking halos before binning works well or not
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
fig_f, (ax_fc, ax_fh)=plt.subplots(2,1,sharex=True) #create two plots with same x-axis for flowline binning
fig_i, (ax_ic, ax_ih)=plt.subplots(2,1,sharex=True) #create two plots with the same x-axis for isochronous binning
flowline_stacked=pd.read_pickle('/home/dbrobins/repos/radfieldcooling/radfieldcooling/saved_dataframes/a0.1451/50_at_mass_10e9_flowline.pkl')
isochronous_stacked=pd.read_pickle('/home/dbrobins/repos/radfieldcooling/radfieldcooling/saved_dataframes/a0.1451/50_at_mass_10e9_isochronous.pkl')

ax_fc.plot(flowline_stacked['temperature'], flowline_stacked['CF_median'], label='all cells median') #Plot the median flowline cooling function
ax_fc.fill_between(flowline_stacked['temperature'], flowline_stacked['CF_25'], flowline_stacked['CF_75'], alpha=0.2, label='all cells spread') #Shade between percentiles of flowline cooling function
ax_fh.plot(flowline_stacked['temperature'], flowline_stacked['HF_median']) #Do the same for flowline heating function
ax_fh.fill_between(flowline_stacked['temperature'], flowline_stacked['HF_25'], flowline_stacked['HF_75'], alpha=0.2)
ax_ic.plot(isochronous_stacked['temperature'], isochronous_stacked['CF_median'], label='all cells median') #Do the same for isochronous binning
ax_ic.fill_between(isochronous_stacked['temperature'], isochronous_stacked['CF_25'], isochronous_stacked['CF_75'], alpha=0.2, label='all cells spread')
ax_ih.plot(isochronous_stacked['temperature'], isochronous_stacked['HF_median'])
ax_ih.fill_between(isochronous_stacked['temperature'], isochronous_stacked['HF_25'], isochronous_stacked['HF_75'], alpha=0.2)

bins=isochronous_stacked['temperature'] #Get temperature bins
flowline_CF_stats=[] #Lists to hold the dataframes for each halo
flowline_HF_stats=[]
isochronous_CF_stats=[]
isochronous_HF_stats=[]
for halo_index in range(50): #Loop through the 20 halos
    flowline_df=pd.read_pickle('/home/dbrobins/repos/radfieldcooling/radfieldcooling/saved_dataframes/a0.1451/mass_10e9_'+str(halo_index)+'_flowline.pkl') #Read in flowline and isochronous dataframes
    flowline_CF=flowline_df['CF_median'].to_numpy()
    flowline_CF_stats.append(flowline_CF.reshape(flowline_CF.shape[0],1))
    flowline_HF=flowline_df['HF_median'].to_numpy()
    flowline_HF_stats.append(flowline_HF.reshape(flowline_HF.shape[0],1))
    isochronous_df=pd.read_pickle('/home/dbrobins/repos/radfieldcooling/radfieldcooling/saved_dataframes/a0.1451/mass_10e9_'+str(halo_index)+'_isochronous.pkl')
    isochronous_CF=isochronous_df['CF_median'].to_numpy()
    isochronous_CF_stats.append(isochronous_CF.reshape(isochronous_CF.shape[0],1))
    isochronous_HF=isochronous_df['HF_median'].to_numpy()
    isochronous_HF_stats.append(isochronous_HF.reshape(isochronous_HF.shape[0],1))
flowline_CF_stats=np.concatenate(flowline_CF_stats, axis=1)
flowline_HF_stats=np.concatenate(flowline_HF_stats, axis=1)
isochronous_CF_stats=np.concatenate(isochronous_CF_stats, axis=1)
isochronous_HF_stats=np.concatenate(isochronous_HF_stats, axis=1)
'''
plt.semilogx(bins, flowline_stacked['CF_median']/np.apply_along_axis(np.nanmedian, 1, flowline_CF_stats), label='Flowline cooling')
plt.semilogx(bins, flowline_stacked['HF_median']/np.apply_along_axis(np.nanmedian, 1, flowline_HF_stats), label='Flowline heating')
plt.semilogx(bins, isochronous_stacked['CF_median']/np.apply_along_axis(np.nanmedian, 1, isochronous_CF_stats), label='Isochronous cooling')
plt.semilogx(bins, isochronous_stacked['HF_median']/np.apply_along_axis(np.nanmedian, 1, isochronous_HF_stats), label='Isochronous heating')
plt.xlabel('Temperature [K]')
plt.ylabel('Ratio (all cells/grouped by halo)')
plt.legend()
plt.savefig('/home/dbrobins/repos/radfieldcooling/radfieldcooling/figures/stacking_ratio_test_50.png')
'''

ax_fc.plot(bins, np.apply_along_axis(np.nanmedian, 1, flowline_CF_stats), linestyle='dashed', label='halo median')
ax_fc.fill_between(bins, np.apply_along_axis(np.nanpercentile, 1, flowline_CF_stats, 25), np.apply_along_axis(np.nanpercentile, 1, flowline_CF_stats, 75), alpha=0.2, label='halo spread')
ax_fh.plot(bins, np.apply_along_axis(np.nanmedian, 1, flowline_HF_stats), linestyle='dashed')
ax_fh.fill_between(bins, np.apply_along_axis(np.nanpercentile, 1, flowline_HF_stats, 25), np.apply_along_axis(np.nanpercentile, 1, flowline_HF_stats, 75), alpha=0.2)
ax_ic.plot(bins, np.apply_along_axis(np.nanmedian, 1, isochronous_CF_stats), linestyle='dashed', label='halo median')
ax_ic.fill_between(bins, np.apply_along_axis(np.nanpercentile, 1, isochronous_CF_stats, 25), np.apply_along_axis(np.nanpercentile, 1, isochronous_CF_stats, 75), alpha=0.2, label='halo spread')
ax_ih.plot(bins, np.apply_along_axis(np.nanmedian, 1, isochronous_HF_stats), linestyle='dashed')
ax_ih.fill_between(bins, np.apply_along_axis(np.nanpercentile, 1, isochronous_HF_stats, 25), np.apply_along_axis(np.nanpercentile, 1, isochronous_HF_stats, 75), alpha=0.2)

ax_fc.set_xscale('log') #Set all axis scales to log
ax_fc.set_yscale('log')
ax_fh.set_xscale('log')
ax_fh.set_yscale('log')
ax_ic.set_xscale('log')
ax_ic.set_yscale('log')
ax_ih.set_xscale('log')
ax_ih.set_yscale('log')
'''
for halo_index in range(20):
    flowline_df=pd.read_pickle('/home/dbrobins/repos/radfieldcooling/radfieldcooling/saved_dataframes/a0.1451/mass_10e9_'+str(halo_index)+'_flowline.pkl')
    isochronous_df=pd.read_pickle('/home/dbrobins/repos/radfieldcooling/radfieldcooling/saved_dataframes/a0.1451/mass_10e9_'+str(halo_index)+'_isochronous.pkl')
    ax_fc.plot(flowline_df['temperature'], flowline_df['CF_median'], linestyle='dotted') #Plot the median flowline cooling and heating functions
    ax_fh.plot(flowline_df['temperature'], flowline_df['HF_median'], linestyle='dotted')
    ax_ic.plot(isochronous_df['temperature'], isochronous_df['CF_median'], linestyle='dotted') #Plot the median isochronous cooling and heating functions
    ax_ih.plot(isochronous_df['temperature'], isochronous_df['HF_median'], linestyle='dotted')
'''

ax_fh.set_xlabel("Temperature T [K]") #Label the axes, create a legend, save the plots                                                                                                     
ax_ih.set_xlabel("Temperature T [K]")
ax_fc.set_ylabel("Flowline $\Lambda$ [erg cm$^3$/s]")
ax_fh.set_ylabel("Flowline $\Gamma$ [erg cm$^3$/s]")
ax_ic.set_ylabel("Isochronous $\Lambda$ [erg cm$^3$/s]")
ax_ih.set_ylabel("Isochronous $\Gamma$ [erg cm$^3$/s]")
ax_fc.legend(loc=4)
ax_ic.legend(loc=2)
fig_f.savefig("/home/dbrobins/repos/radfieldcooling/radfieldcooling/figures/two_stacking_methods_flowline_50.png")                                                                            
fig_i.savefig('/home/dbrobins/repos/radfieldcooling/radfieldcooling/figures/two_stacking_methods_isochronous_50.png')

'''
fig_f.savefig("/home/dbrobins/repos/radfieldcooling/radfieldcooling/figures/test_stacking_flowline.png")
fig_i.savefig('/home/dbrobins/repos/radfieldcooling/radfieldcooling/figures/test_stacking_isochronous.png')
'''
