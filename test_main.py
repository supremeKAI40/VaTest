import Nearest_Neighbour as test

#testing for toi_1973 sector 38
tics=[47617161]
sector=[38]
transitlist=[] #Should be transit times in the form of list !!NOT TESTED YET

'''path of where all sector wise tics with RA and dec are situated
!!!Careful on path to be Windows/Linux format according to user's OS
it is situated in tess mit site online @https://tess.mit.edu/observations/target-lists/
'''
indir="C:/Personal/Internships/ExoplanetStatisticalValidation/lattetest/test_automation"


# --------------------------------------------
# -----------------testing ------------------
# --------------------------------------------

'Calling tic id list function'
ticids, distance, target_ra, target_dec=test.nn_ticids(indir,sector,tics)
print("Neighbour tics: ",ticids,"\n","corresponding distance: ",distance,"\n","RA: ",target_ra," & ","Dec: ",target_dec)

'Calling tic download function'
alltime, allflux, all_md, alltimebinned, allfluxbinned, outtics, tessmag_list, distance = test.download_data_neighbours(indir, sector, ticids, distance, binfac = 5)


'Calling tic neighbour plot function'
test.plot_nn(tics,indir,alltime,allflux,alltimebinned, allfluxbinned,transitlist,ticids,tessmag_list,distance,True)