import matplotlib.pyplot as plt
import numpy as np
import os

def plot_centroids(lc, transit_time=None, window=None, i = None):
        
    # bin the data 
    lc_bin = lc.bin(7/60/24)
    
    # generate a mask so that we only see the times around the transit event
    if (window == None) or (transit_time == None):
        transit_mask = lc_bin.time.value > 0
    else:
        transit_mask = (lc_bin.time.value > transit_time - window) & (lc_bin.time.value < transit_time + window)
        
    # make a plot with three panels so that we can see the lightcurve and the centroid positions 
    fig, ax = plt.subplots(3,1, figsize = (8,5), sharex= True)
    
    # plot the lightcurve in the top panel (in orange) 
    ax[0].plot(lc_bin.time.value[transit_mask], lc_bin.sap_flux.value[transit_mask], color = 'darkorange', lw = 0, marker = '.', ms = 3)
    
    # plot the centroid motions in the column direction in the middle panel
    ax[1].plot(lc_bin.time.value[transit_mask], lc_bin.mom_centr1.value[transit_mask]-np.nanmean(lc_bin.mom_centr1.value[transit_mask]), color = 'black', lw = 0, marker = '.', ms = 2, alpha = 0.8)
    ax[1].plot(lc_bin.time.value[transit_mask], lc_bin.pos_corr1.value[transit_mask]-np.nanmean(lc_bin.pos_corr1.value[transit_mask]), color = 'red', lw = 0, marker = '.', ms = 2, alpha = 0.8)
    
    # plot the centroid motions in the row direction in the middle panel
    ax[2].plot(lc_bin.time.value[transit_mask], lc_bin.mom_centr2.value[transit_mask] - np.nanmedian(lc_bin.mom_centr2.value[transit_mask]), color = 'black', lw = 0, marker = '.', ms = 2, alpha = 0.8, label = 'Brightness motion')
    ax[2].plot(lc_bin.time.value[transit_mask], lc_bin.pos_corr2.value[transit_mask]- np.nanmedian(lc_bin.pos_corr2.value[transit_mask]), color = 'red', lw = 0, marker = '.', ms = 2, alpha = 0.8, label = 'Satellite motion')
    
    if (transit_time != None):
        # draw a vertical line at the time of the transit event 
        ax[0].axvline(transit_time, color = 'grey', zorder = -1)
        ax[1].axvline(transit_time, color = 'grey', zorder = -1)
        ax[2].axvline(transit_time, color = 'grey', zorder = -1)
    
    # label the axes
    ax[0].set_ylabel("Flux")
    ax[1].set_ylabel("Column position")
    ax[2].set_ylabel("Row position")
    plt.xlabel("Time (TJD)")
    
    if (window != None) and (transit_time != None):   
        plt.xlim(transit_time - window, transit_time + window)
        
    plt.tight_layout()
    plt.legend()
    plt.show()
    
    #fig.savefig(os.path.join('', 'TOI_181_C02'+str(i)+'.jpeg'), format='jpeg', dpi=600, bbox_inches = 'tight')