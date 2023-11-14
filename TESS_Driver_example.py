#TESS EBLM Driverimport sys
import sys
sys.path.append("PackageFiles")
from ESO_harps_download3 import ESOdownload
from HARPSrvPipline3 import RVharps
from GeneratingFreshProFile import MakeProFile
from MakingFreshPriorFile import MakeNewPriorFile, UpdatePriorFile
from TESS_pipeline_wotan_alisonchanges import startTess
from MakeSbatchFile import MakeSbatch
from lightcurve_pipeline import create_light_curve
import os
import pandas as pd
import numpy as np
import lightkurve as lk
import time


#I THINK I SHOULD JUST HARD CODE MAKING A NEW START OF THE FILE
#command =f"cp PackageFiles/makepriors.sh makepriors.sh"
#os.system(command)

with open("errorlogEBLM.txt", 'w') as f:
		f.writelines("TESS Pipeline Error Log \n")



def start(shortname, ticID, maxi = 0 , starval=0):
	'''
	Shortname: The common name of your system
	ticID: The TIC ID
	maxi: the maximum number of systems you want set up files for. The default is 1.
	starval: If looping over several systems, the index of the current system The default is index zero.


	'''
	print(shortname)	
	command =f"pwd"
	os.system(command)

	#Checking if transit data or target information is available.

	try:
		ticName = "TIC "+str(ticID)
		lc = lk.search_lightcurve(ticName,author='SPOC',mission="TESS",exptime=120)
		short_mission = [m[-2:] for m in lc.mission]
	
	except:
		print(f"{starval} - {shortname}: No TESS Observations")
		with open("errorlogEBLM.txt", 'a', encoding='utf-8') as f:
			f.writelines(f"{starval} - {shortname}: No TESS Observations \n")
			return
	
	try:
		url = f'https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/{ticID}/info/?tce=1&sector={short_mission[0]}'
		params = pd.read_json(url,orient='index')
		duration_1sig = (params['TDUR'][1] + params['INDUR'][1])/24 

	except:
		print(f"{starval} - {shortname}: No TESS Observations")
		with open("errorlogEBLM.txt", 'a', encoding='utf-8') as f:
			f.writelines(f"{starval} - {shortname}: {url} \n")
			f.writelines(f"{starval} - {shortname}: No information listed in TIC \n")
			return


	#Setting up file pathways
	startdir = os.getcwd()
	if os.path.exists('./'+shortname):
		os.chdir('./'+shortname)
	else:
		os.mkdir('./'+shortname)
		os.chdir('./'+shortname)




	tic = time.perf_counter()

	#print("Downloading data from the ESO")
	#ESOdownload(ticID,shortname,'HARPS')
	#ESOdownload(ticID,shortname,'ESPRESSO')

	#print("Cleaning RV data")
	#RVharps(ticID,shortname,'HARPS')
	#RVharps(ticID,shortname,'ESPRESSO')

	print("Downloading and Flattening TESS Observations")
	#startTess(ticID,shortname)
	for m in short_mission:
		create_light_curve(f'TIC {ticID}', 'SPOC', int(short_mission[0]), auto = True, exposure=120, targetname = shortname, 
			save = True, plot = False, planet='b', published = False, qualityFlag = True)

	print("Making the Pro File")
	MakeProFile(ticID,shortname)

	print("Making the Prior File Commands")
	MakeNewPriorFile(ticID,shortname,starval,maxi)

	print("Making the Sbatch Files")
	MakeSbatch(shortname)



	toc = time.perf_counter()
	print(f"Completed in {toc - tic:0.4f} seconds")

	#Run this section after you have the exofast generated prior files
	#print("Updating the Prior File")
	#UpdatePriorFile(shortname,ticID) 


	os.chdir('..')
	print()
	#command =f"pwd"
	#os.system(command)



shortname = "HD 209458b"
shortname = shortname.replace(" ","")
ticID = 420814525

print()
print()
print()
start(shortname,ticID)


	

