import requests
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.io.fits as pf
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np


def nn_ticids(indir, transit_sec, tic):
    '''
    function to find the TIC IDs of the 6 nearest neighbour stars that were observed by TESS (short cadence).
    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    transit_sec  : list
        list of the sectors in which the peaks appear.
    tic : str
        TIC (Tess Input Catalog) ID of the target
    Returns
    -------
    ticids  : list
        TIC IDs of the 6 TESS target pixels file stars that are closest to the target.
    target_ra  : float
        Right Ascension of the target
    target_dec  : float
        Declination of the target
    '''

    neighbours_sector = transit_sec[0]
'''
    if neighbours_sector < 10:
        download_sector = "00{}".format(neighbours_sector)
    else:
        download_sector = "0{}".format(neighbours_sector)

    # load the data file of the sector that the first marked transit appears in
    # sort the list to be ordered by camera, RA and then Dec
    #tic_list = pd.read_csv("{}/data/all_targets_S{}_v1.txt".format(indir,download_sector)).sort_values(['Camera', 'RA', 'Dec']).reset_index()
'''    
    tic_list_all = pd.read_csv("{}/data/all_targets_list.txt".format(indir)).sort_values(['sec','Camera', 'RA', 'Dec']).reset_index()
    
    tic_list=tic_list_all[(tic_list_all['sec']==neighbours_sector[0])].reset_index()


    # acess the information for the target stars
    target = tic_list.loc[tic_list['TICID'] == float(tic[0])]
    tic_idx = target.index[0]  # get the index of the target in the list

    # get the RA and Dec of the target star - this will be returned but the function
    target_ra = float(target['RA'])
    target_dec = float(target['Dec'])

    # ------------

    # make a list of the closest stars to the target (only stars that are TESS target stars are considered)
    # ensure that the TIC ID is not at the end or begining of the tic list - otherwise we can't take 100 from wither side.
    if tic_idx < 101:
        tic_list_close = tic_list[0:tic_idx + 101]
    elif tic_idx > (len(tic_list) + 1):
        tic_list_close = tic_list[tic_idx - 100 :len(tic_list)]
    else:
        tic_list_close = tic_list[tic_idx - 100:tic_idx + 101]
    # ------------

    # function to alculated the angular separation from the target to the stars to find the nearest neighbours
    def star_sep(row, ra, dec):

        ra2 = float(row['RA'])
        dec2 = float(row['Dec'])

        c1 = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
        c2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')
        sep = c1.separation(c2)

        return sep.arcminute

    # make a new column in the pandas dataframe of the separation
    tic_list_close['dist'] = tic_list_close.apply(star_sep, args = (target_ra,target_dec), axis=1)

    # Select the 6 nearest neightbours and creare a list of their tic IDs and distance to them.
    closest_tic = tic_list_close.sort_values('dist')[0:6]
    ticids = closest_tic['TICID'].tolist()
    distance = closest_tic['dist'].tolist()

    return ticids, distance, target_ra, target_dec

#Function to rebin the nearest neighbours
def rebin(arr,new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)

def download_data_neighbours(indir, sector, tics, distance, binfac = 5):
    ''''
    Dowloand the data for the 6 nearest neightbour target pixel LC - in the future want to
    take this so it used the FFI data as this will be closer targets and a better diagnostic.
    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (defaul = "./LATTE_output")
    sector  :  list or str
        list of the sectors that we want to analyse. If 'all', all the sectors in whic the target appears will be downloaded.
    distance:  list
        list of the distance to the nearest neighbouring tics from the target (in arcminutes)
    tics : str
        TIC (Tess Input Catalog) ID of the target
    binfac  :  int
        The factor by which the data should be binned. Default = 5 (which is what is shown on PHT)
    Returns
    ------
    alltime  :  list
        times (not binned)
    allflux  :  list
        normalized flux (not binned)
    all_md  :  list
        times of the momentum dumps
    alltimebinned  :  list
        binned time
    allfluxbinned  :  list
        normalized binned flux
    transit_list  :  list
        list of all the marked transit-events
    outtics  :  list
        the tic ids of the nearby targets
    tessmag_list  :  list
        the magnitudes of the nearby targets
    distance  :  list
        the distances to the nearby targets
    '''

    sb.set(style='ticks')
    sb.set_color_codes('deep')

    #try:
    #    lc_sec = np.genfromtxt('{}/data/tesscurl_sector_{}_lc.sh'.format(indir, str(sector)), dtype = str)
    #except:
    #    print ("Sector {}  has not yet been osberved - come back to this target later.".format(sector))

    alltime_neighbours = []
    allflux_neighbours = []
    all_md_neighbours = []
    alltimebinned_neighbours = []
    allfluxbinned_neighbours = []
    outtics = []

    dwload_link = []

    sector_int=int(sector[0])

    infile = pd.read_csv("{}/data/sector_download_codes.txt".format(indir), delimiter = ' ', names = ['sec', 'first', 'second'], comment = '#')
    for tic in tics:
        
        this_sector_code = infile[infile.sec == sector_int]

        sector_str = str(sector[0])

        download_url = (
            "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess"
            + str(this_sector_code['first'].values[0]).rjust(13, "0")
            + "-s"
            + sector_str.rjust(4, "0")
            + "-"
            + str(tic).rjust(16, "0")
            + "-"
            + str(this_sector_code['second'].values[0]).rjust(4, "0")
            + "-s_lc.fits")

        dwload_link.append(download_url)


    alltimebinned = []
    allfluxbinned = []

    start_sec = []
    end_sec = []

    alltime = []
    allflux = []
    all_md = []

    tessmag_list = []

    for num,lcfile in enumerate(dwload_link):

        print ("Downloading nearest neighbour   {}   of   {}....".format(num + 1, len(dwload_link), end ='' ))
        try:
            response = requests.get(lcfile)

            # open the file using the response url
            lchdu  = pf.open(response.url) # this needs to be a URL - not a file
            outtics.append(int(lchdu[0].header['TICID']))

            #Open and view columns in lightcurve extension
            lcdata = lchdu[1].data
            lchdu[1].columns

            f02 = lcdata['PDCSAP_FLUX']
            quality = lcdata['QUALITY']
            time    = lcdata['TIME']
            f0     = lcdata['SAP_FLUX']
            fbkg     = lcdata['SAP_BKG']

            med = np.nanmedian(f02)
            f1 = f02/med

            x1      = lcdata['MOM_CENTR1']
            x1      -= np.nanmedian(x1)
            y1      = lcdata['MOM_CENTR2']
            y1      -= np.nanmedian(y1)
            x2      = lcdata['POS_CORR1']
            x2      -= np.nanmedian(x2)
            y2      = lcdata['POS_CORR2']
            y2      -= np.nanmedian(y2)
            l       = (quality>0)
            l2      = (quality<=0)

            sec     = int(lchdu[0].header['SECTOR'])

            flux     = lcdata['SAP_FLUX']

            tessmag = lchdu[0].header['TESSMAG']

            lchdu.close()

            # binned data
            N       = len(time)
            n       = int(np.floor(N/binfac)*binfac)
            X       = np.zeros((2,n))
            X[0,:]  = time[:n]
            X[1,:]  = f1[:n]
            Xb      = rebin(X, (2,int(n/binfac)))

            time_binned    = Xb[0]
            flux_binned    = Xb[1]

            mom_dump = np.bitwise_and(quality, 2**5) >= 1

            alltime.append(list(time)) #[~bad_data]
            allflux.append(list(f1)) #[~bad_data])
            all_md.append(list(time[mom_dump]))

            alltimebinned.append(list(time_binned))
            allfluxbinned.append(list(flux_binned))

            start_sec.append([time[0]])
            end_sec.append([time[-1]])
            tessmag_list.append(tessmag)

            print("done.")
        except:
            continue
        
    return alltime, allflux, all_md, alltimebinned, allfluxbinned, outtics, tessmag_list, distance


# plot nearest neighbour LCs
def plot_nn(tic, indir, alltime_nn, allflux_nn, alltimebinned_nn, allfluxbinned_nn, transit_list, outtics, tessmag_list, distance, args):

    '''
    Plot the lighcurves of the 6 nearest neighbours to the target.
    Distance to the targets is calculated using the RA and DEC of the targte and adding them in quadrature
    Parameters
    ----------
    tic : str
        TIC (Tess Input Catalog) ID of the target
    indir   :   str
        path to where the data will be saved (default = "./LATTE_output")
    alltime_nn  :  list
        times for all the nearest neighbours (not binned)
    allflux_nn  :  list
        normalized flux times for all the nearest neighbours (not binned)
    alltimebinned_nn  :  list
        binned time for all the nearest neighbours
    allfluxbinned_nn  :  list
        normalized binned flux for all the nearest neighbours
    transit_list  :  list
        list of all the marked transits
    outtics  :  list
        the tic IDs of the 6 nearest neighbours
    Returns
    -------
        Plot of the 6 nearest neighbours. The vertical line indicated the location(s) of the marked transit.
    '''

    fig, ax = plt.subplots(len(alltime_nn), 1, figsize=(13,8), sharex=True, gridspec_kw={'hspace': 0})
    plt.tight_layout()

    colors = ['r', 'darkorange', 'gold', 'seagreen', 'royalblue', 'navy','magenta' ]
    colors2 = ['k', 'k', 'k', 'k', 'k', 'grey','k' ]
    for i in range(0,(len(alltime_nn))):

        for line in (transit_list):
            ax[i].axvline(line, color = 'k', linewidth = 2.2, alpha = 1, linestyle = '-')
        if str(outtics[i]) == str(tic):
            ax[i].plot(alltime_nn[i], np.array(allflux_nn[i]), color = colors[i], label = "*** {}  Tmag = {:3f}***".format(tic, tessmag_list[i]), marker = '.', ms = 2, linewidth = 0)

        else:
            ax[i].plot(alltime_nn[i], np.array(allflux_nn[i]), color = colors[i], label = "{}  Tmag = {:3f}   d = {:3f} arcmins".format(outtics[i], tessmag_list[i], distance[i]), marker = '.', ms = 2, linewidth = 0)

        ax[i].plot(alltimebinned_nn[i], np.array(allfluxbinned_nn[i]), color = colors2[i], marker = '.', ms = 1, linewidth = 0)

        ax[i].legend(loc = 1)

    ax[0].set_title("LCs of Nearby Stars", fontsize = 13)
    ax[0].set_title("LCs of Nearby Stars", fontsize = 13)
    ax[len(alltime_nn) - 1].set_xlabel("Time (BJD-2457000)", fontsize = 13)
    ax[int(len(alltime_nn)/2)].set_ylabel("Normalised Flux", fontsize = 13)
    ax[0].set_xlim(np.nanmin(alltime_nn), np.nanmax(alltime_nn))

    plt.xlim(np.nanmin(alltime_nn), np.nanmax(alltime_nn))

    plt.savefig('{}/Output/{}_nearest_neighbours.png'.format(indir, tic[0]), format='png', bbox_inches='tight')


# plot the image cut-out of the nearby stars
def plot_cutout(image):
    """
    Plot image cut out of the target.
    """
    plt.imshow(image, origin = 'lower', cmap = plt.cm.YlGnBu_r,
           vmax = np.percentile(image, 92),
           vmin = np.percentile(image, 5))

    plt.grid(axis = 'both',color = 'white', ls = 'solid')