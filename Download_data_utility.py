import os
import shutil
import requests
import numpy as np
import pandas as pd

def tp_files(indir):
    '''
    Function to download all of the TIC list and RA/Dec data that we want to the local computer.
    These are needed for the nearest neighbour analysis.
    Parameters
    ----------
    indir   :   str
        path to where the data will be saved (default = "")
    '''
    #create the folder where nearest neighbour analysis file will look for data
    if not os.path.exists("{}/data".format(indir)):
        os.mkdir("{}/data".format(indir))
        
    if not os.path.exists("{}/data/all_targets_list.txt".format(indir)):
        
        empty_dataframe = pd.DataFrame(columns=['TICID', 'Camera', 'CCD', 'Tmag', 'RA', 'Dec', 'sec'])

        empty_dataframe.to_csv("{}/data/all_targets_list.txt".format(indir), index=False)

        first_sec = 0 # start with sector 1 but this has to be 0 because the next step of the code adds one (needs to be like this otherwise it will dowload the last sector multiple times when re run)
        print ("We also need to have a record of all of the TIC ID available in each Sectors \nStarting with Sector 1")

    else:
        all_sectors_file = pd.read_csv("{}/data/all_targets_list.txt".format(indir), delimiter = ',')
        exist = all_sectors_file['sec'].max()
        first_sec = (np.max(exist))
        
    for sec in range(first_sec+1,500):  # max sector is varying for tess

        if sec < 10:
            download_sector = "00{}".format(sec)
        else:
            download_sector = "0{}".format(sec)

        target_list = "https://tess.mit.edu/wp-content/uploads/all_targets_S{}_v1.txt".format(download_sector)

        r_target_list = requests.get(target_list) # create HTTP response object

        if r_target_list.status_code == 404:
            print ("Target lists only available up to Sector {} -- try downloading more data later".format(sec))
            break
        
        with open("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector), 'wb') as f:
            f.write(r_target_list.content)

        #save the response as csv add sector as column and then add to all target file
        this_sector_file = pd.read_csv("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector), comment = '#', delimiter = '\t', names = ['TICID', 'Camera', 'CCD', 'Tmag', 'RA', 'Dec'])
        this_sector_file['sec'] = sec
        this_sector_file.to_csv("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector), index=False)

        # also open the latest version of the 'all' file so that we can add the new part to it... again not the best way
        all_sectors_file = pd.read_csv("{}/data/all_targets_list.txt".format(indir, download_sector), comment = '#', delimiter = ',')

        # merge the dataframes
        full_data_frame_allsectors = pd.concat([all_sectors_file, this_sector_file])

        # save it again
        full_data_frame_allsectors.to_csv("{}/data/all_targets_list.txt".format(indir, download_sector), index=False)
        
        print("Finished adding target list text file for sector {}".format(sec))
        if os.path.exists("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector)):
            os.remove("{}/data/all_targets_S{}_v1.txt".format(indir, download_sector))        
         
        ##Adding sector download codes
        LC_url = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_lc.sh".format(sec)
        r_LC = requests.get(LC_url) # create HTTP response object

        if r_LC.status_code == 404:
            print ("You're all caught up with your codes up to Sector {} -- try downloading more data later".format(sec-1))
            break


        with open("{}/data/sector_download_codes.txt".format(indir), 'a') as f:
                
                #open the first line of the file in order to get the codes that are requred to download data from that sector
                
                code1 = (str(r_LC.content[0:200]).split('-')[4][6:])
                code2 = (str(r_LC.content[0:200]).split('-')[7])

                string = "{} {} {}\n".format(str(sec), code1, code2 )

                # write the contents of the response (r.content)
                # to a new file in binary mode.

                f.write(string)
                print("Finished adding sector codes for sector {}".format(sec))