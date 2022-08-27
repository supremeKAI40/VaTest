import lightkurve as lk  # LightKurve Python Package to View & Process the Transit Light Curves
import numpy as np       # Numpy for Scientific Computing
import matplotlib.pyplot as plt # Matplotlib for Plotting the Garphs
import os
import pandas as pd

##inputs sample to refer during build

'''
input list
    toi: as array i.e [181]
    Sector: as array i.e [38]
    Indir: the file path used to access sector wise download code
'''

def combine_lc(tois,sectors,indir):
    ##to download the sector wise lightkurve
    search_result = lk.search_lightcurve('TOI '+str(tois[0]), author='SPOC',exptime='short')
    lc_collection=search_result.download_all()
    Combined_TOI=lc_collection.stitch()
    print(search_result)

    pg = Combined_TOI.to_periodogram(method = 'bls', period = np.arange(1,20,0.01))
    pg.plot()

    time_period = pg.period_at_max_power
    total_transit_time = pg.transit_time_at_max_power
    transit_duration = pg.duration_at_max_power
    print("Parameters from first pass \n")
    print(" Orbital Period (d) : ", time_period, "\n Epoch Time (BTJD) : ", total_transit_time,\
        "\n Transit Duration (d) : ", transit_duration)

    ##masking to detrend and get more accurate transit avoiding stellar activity

    masked_lc = Combined_TOI.create_transit_mask(period = time_period, 
                                             transit_time = total_transit_time , 
                                             duration = transit_duration.value+0.05)

    # Detrending the out of transit part
    flat_lc = Combined_TOI.flatten(mask = masked_lc)

    # Create the Periodogram of Detrended Light Curve
    pg_flatten = flat_lc.to_periodogram(method = 'bls', period = np.arange(1,20,0.001))
    pg_flatten.plot()
    plt.show()
    # Print the evaluated values
    time_period = pg_flatten.period_at_max_power
    transit_time = pg_flatten.transit_time_at_max_power
    transit_duration = pg_flatten.duration_at_max_power
    print("\nParameters from detrended periodogram\n")
    print(" Orbital Period (d) : ",time_period,"\n","Epoch Time (BTJD) : ",transit_time,\
      "\n","Transit Duration (d) : ", transit_duration)

toi= [181]
sector=[2,29]
indir="C:/Personal/Internships/ExoplanetStatisticalValidation/VATEST_package/"

combine_lc(toi,sector,indir)