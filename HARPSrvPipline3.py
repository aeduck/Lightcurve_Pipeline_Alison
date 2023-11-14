import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from astropy.time import Time
from urllib.request import urlopen
import tarfile
import os
import pandas as pd
from io import StringIO
import sys
import datetime


def RVharps(ticID,shortname,inst):
    allRV1_values = []
    allRV2_values = []
    allespresso_values = []

    if os.path.exists('./archive'):
        os.chdir('./archive')
    else: 
        print("No archive!")
        return
    
    #print(os.listdir())
    for x in os.listdir():
        neededfile = False
        print(x)
        if inst=='HARPS' and x.endswith(".tar"):
            neededfile = True
            tarfilename=x
            print(tarfilename)
            rv_table = []
            error = False

            file_obj= tarfile.open(tarfilename,"r")
            namelist=file_obj.getnames()

            for name in namelist:
                    if "ccf" in name and name.endswith("_A.fits"):
                        print(name)
                        rv_table = name
            file=file_obj.extractfile(rv_table)

        elif inst=='ESPRESSO' and x.endswith(".fits"):
            neededfile = True
            file = x
        else:
            continue
     

        
        hdul = fits.open(file)
        if inst == 'HARPS':
            bjdtime = hdul[0].header['HIERARCH ESO DRS BJD']
            drs_ccf_rvc = hdul[0].header['HIERARCH ESO DRS CCF RVC']*1000
            try:
                drs_ccf_noise = hdul[0].header['HIERARCH ESO DRS CCF NOISE']*1000
                print(hdul[0].header['HIERARCH ESO DRS CCF NOISE'])
                print(drs_ccf_noise)
            except:
                print(f"Data are from{hdul[0].header['HIERARCH ESO TPL START']}")
                print(f"Skipping observations")
                
        elif inst == 'ESPRESSO':
            bjdtime = hdul[0].header['HIERARCH ESO QC BJD']
            drs_ccf_rvc = hdul[0].header['HIERARCH ESO QC CCF RV']*1000
            drs_ccf_noise = hdul[0].header['HIERARCH ESO QC CCF RV ERROR']*1000
            print(drs_ccf_noise)
                
        
        if neededfile == True and inst == 'HARPS':
            year = int(x[4:8])
            month = int(x[9:11])
            day = int(x[12:14])

            d1 = datetime.datetime(year, month, day)
            d2 = datetime.datetime(2015, 6, 15)


            if d1 < d2:
                allRV1_values.append([bjdtime,drs_ccf_rvc,drs_ccf_noise])
            else:
                allRV2_values.append([bjdtime,drs_ccf_rvc,drs_ccf_noise])
        
        if neededfile and inst == 'ESPRESSO':
            allespresso_values.append([bjdtime,drs_ccf_rvc,drs_ccf_noise])

        print()
        



    if not os.path.exists('../Data'):
        os.mkdir('../Data')


    if allRV1_values:
        np.savetxt(f'../Data/{shortname}.HARPS.HARPS1.rv', allRV1_values, delimiter =" ",fmt ='% s')
    if allRV2_values:
        np.savetxt(f'../Data/{shortname}.HARPS.HARPS2.rv', allRV2_values, delimiter =" ",fmt ='% s')
    if allespresso_values:
        np.savetxt(f'../Data/{shortname}.ESPRESSO.ESPRESSO.rv', allespresso_values, delimiter =" ",fmt ='% s')

    os.chdir('..')




if __name__ == "__main__":

    shortname = "HD209458"
    ticID = 420814525
    
    inst = 'ESPRESSO'
    RVharps(ticID,shortname,inst)
    inst = 'HARPS'
    RVharps(ticID,shortname,inst)

