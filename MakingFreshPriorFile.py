
import numpy as np
import pandas as pd
import sys
import os
import math
import lightkurve as lk
from astropy.constants import R_jup, R_sun


def MakeNewPriorFile(ticID,shortname,starval=0,maxi=0,pathexo = "/home/duck.18/idl/EXOFASTv2/Benchmark",email=""):

	"""
  	Creates a new prior file for EXOFASTv2.

  	Args:
    	starval: The current star number.
    	maxi: The maximum number of stars in your list.
    	shortname: The short name of the system.
    	ticID: The TIC ID of the system.
	    pathexo: string of the path where you want all your files to live (string)
		email: email adress for Unity to send emails to

  	Returns:
    	None.
  	"""
	  
	#print(shortname)
	#print(os.system("pwd"))

	ticName = "TIC "+str(ticID)

	#idlcommand = "mkticksed,'382188882',priorfile="J0525.priors",sed="J0525.sed" "

	idlcommand = f"mkticsed, \'{ticID}\',priorfile=\'{shortname}.priors\',sed=\'{shortname}.sed\'"

	print(idlcommand)

	boilerplate =["#!/usr/bin/env bash\n","#SBATCH --time=00:30:00 \n","#SBATCH --nodes=1 \n","#SBATCH --ntasks-per-node=11\n","#SBATCH --mem=16g\n",f"#SBATCH --job-name=priors\n","#SBATCH --mail-type=ALL\n","#SBATCH --mail-user={email}\n", "source ~/source.sh\n", "ml idl"]
	
	
	if starval == 0:
		with open(f'../makepriors2.sh', 'w', encoding='utf-8') as file:
			file.writelines(boilerplate)
			file.writelines(f"\n")
			file.writelines(f"cd {pathexo} \n")
			file.writelines(f"\n")
		file.close()
	with open(f'../makepriors2.sh', 'a', encoding='utf-8') as file:
		file.writelines(f'idl -e "{idlcommand}" \n')
		file.writelines(f'mv {shortname}.priors {shortname}/SetupFiles/{shortname}.priors \n')
		file.writelines(f'mv {shortname}.sed {shortname}/Data/{shortname}.sed \n')
	file.close()

	if starval == maxi:
		with open(f'../makepriors2.sh', 'a', encoding='utf-8') as file:
			file.writelines("wait")
			file.writelines('\n')
		file.close()
	print(os.path.realpath(file.name))
	print()




	

def UpdatePriorFile(shortname,ticID):
	"""
	Updates the prior file for EXOFASTv2 using the values from the TESS Input Catalog, Gaia EDR3, and the Exoplanet Archive.

	Args:
		shortname: The short name of the system.
		ticID: The TIC ID of the system.

	Returns:
		None.
  	"""
	GaiaDR3 = True


	#url = r"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,tic_id,default_flag,pl_orbper,pl_ratdor,pl_ratror,pl_trandep,pl_radj,st_rad,pl_orbsmax,st_teff,st_tefferr1,st_teffer2,st_met,st_meterr1,st_meterr2+from+ps+where+tic_id=%27TIC%20"+str(ticID)+"%27+and+default_flag=1+&format=csv"
	url =r"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,tic_id,default_flag,pl_orbper,pl_ratdor,pl_ratror,pl_trandep,pl_radj,st_rad,pl_orbsmax,st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2+from+ps+where+tic_id=%27TIC%20"+str(ticID)+"%27+and+default_flag=1+&format=csv"
	print(url)
	params = pd.read_table(url,sep=',')

	oneparams = params.iloc[0]

	print(oneparams)

	print()

	print(f"Name: {oneparams['pl_name']}")
	print(f"Period in Days: {oneparams['pl_orbper']}")
	print(f"R_planet/R_star: {oneparams['pl_ratror']}")
	print(f"Transit Depth [%]: {oneparams['pl_trandep']}")
	print(f"pl_radj/st_rad: {(oneparams['pl_radj']*R_jup)/(oneparams['st_rad']*R_sun): .3f}")
	print(f"a/R_star: {oneparams['pl_ratdor']}")
	print(f"a [au]: {oneparams['pl_orbsmax']}")
	print(f"Teff [k]: {oneparams['st_teff']}")
	print(f"Teff err1[k]: {oneparams['st_tefferr1']}")
	print(f"Teff err2[k]: {oneparams['st_tefferr2']}")
	print(f"Metalliticy [dex]: {oneparams['st_met']}")
	print(f"Metalliticy err1[dex]: {oneparams['st_meterr1']}")
	print(f"Metalliticy err2[dex]: {oneparams['st_meterr2']}")



	teff_eblm = oneparams['st_teff']
	teff_eblm_err = (np.abs(oneparams['st_tefferr1']) + np.abs(oneparams['st_tefferr2']))/2
	if np.isnan(teff_eblm_err):
		teff_eblm_err=''
	feh_eblm = oneparams['st_met']
	if np.isnan(float(feh_eblm)):
		feh_eblm=''
	feh_eblm_err = (np.abs(float(oneparams['st_meterr1'])) + np.abs(float(oneparams['st_meterr2'])))/2
	if np.isnan(feh_eblm_err):
		feh_eblm_err=''

	
	command =f"pwd"
	os.system(command)

	#with open(f'./{shortname}/SetupFiles/{shortname}.priors', 'r') as file:
	with open(f'./SetupFiles/{shortname}.priors', 'r') as file:
	#with open(f'{shortname}.priors', 'r') as file:

		data = file.readlines()
		noteff = True
		nofeh = True


		for i, line in enumerate(data):

			if 'teff' in line:
				print("Found Teff")
				data[i] = f"teff {teff_eblm} {teff_eblm_err} \n"
				noteff = False

			elif 'feh' in line:
				print("Found feh")
				if feh_eblm == '':
					data[i] = f"\n"
				else:
					data[i] = f"feh {feh_eblm} {feh_eblm_err} \n"
				nofeh = False

			if not noteff and not nofeh:
				break




	if noteff or nofeh:	
		if noteff:
			print("Did not find Teff")
			for i, line in enumerate(data):
				if '##############' in line:
					linenum = i
					break
			data.insert(linenum-1, f"teff {teff_eblm} {teff_eblm_err} \n")

		if nofeh:
			print("Did not find feh")
			for i, line in enumerate(data):
					if '##############' in line:
						linenum = i
						break
			if feh_eblm == '':
				data.insert(i-1, f"feh 0.0 1.0 \n")
			else:
				data.insert(i-1, f"feh {feh_eblm} {feh_eblm_err} \n")






	if GaiaDR3:
		parallaxcoutner = 0
		for i,line in enumerate(data):
			if '#parallax' in line:
				data[i] = data[i].replace("#parallax","parallax",1)

			elif 'parallax' in line:
				if parallaxcoutner == 1:
					data[i] = data[i].replace("parallax","#parallax",1)
				parallaxcoutner = parallaxcoutner +1


	#with open(f'./{shortname}/SetupFiles/{shortname}.priors.1', 'w', encoding='utf-8') as file:
	with open(f'./SetupFiles/{shortname}.priors.1', 'w', encoding='utf-8') as file:
	#with open(f'{shortname}.priors.1', 'w', encoding='utf-8') as file:
		file.writelines(data)


	if GaiaDR3:
		#with open(f'./{shortname}/SetupFiles/{shortname}.sed', 'r') as file:
		with open(f'./Data/{shortname}.sed', 'r') as file:
		#with open(f'{shortname}.sed', 'r') as file:
			data = file.readlines()
			for i,line in enumerate(data):
				if '# Gaia_G_EDR3' in line:
					data[i] = data[i].replace("# Gaia_G_EDR3","Gaia_G_EDR3")
				if '#Gaia_BP_EDR3' in line:
					data[i] = data[i].replace("#Gaia_BP_EDR3","Gaia_BP_EDR3")
				if '#Gaia_RP_EDR3' in line:
					data[i] = data[i].replace("#Gaia_RP_EDR3","Gaia_RP_EDR3")

				if 'GaiaBP ' in line:
					print("found GaiaBP ")
					data[i] = data[i].replace("GaiaBP ","#GaiaBP ")
				if 'Gaia ' in line:
					print("found Gaia ")
					data[i] = data[i].replace("Gaia ","#Gaia ")
				if 'GaiaRP ' in line:
					print("found GaiaRP ")
					data[i] = data[i].replace("GaiaRP ","#GaiaRP ")


		#with open(f'./{shortname}/Data/{shortname}.sed', 'w', encoding='utf-8') as file:
		with open(f'./Data/{shortname}.sed', 'w', encoding='utf-8') as file:
		#with open(f'{shortname}.sed', 'w', encoding='utf-8') as file:
			file.writelines(data)


	

	ticName = "TIC "+str(ticID)
	lc = lk.search_lightcurve(ticName,author='SPOC',mission="TESS",exptime=120)
	short_mission = [m[-2:] for m in lc.mission]
	url = f'https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/{ticID}/info/?tce=1&sector={short_mission[0]}'
	print(url)
	params = pd.read_json(url,orient='index')


	periodEBLM = float(params["TPERIOD"][1])
	tc_dur = float(params["TDUR"][1])/24
	tc_tic = float(params["TEPOCH"][1])+float(params["BJDREFI"][1])
	RadRatio = float(params["RADRATIO"][1])
	RadBool = True
	print(f"{periodEBLM}\n {tc_dur}\n {tc_tic} \n{RadRatio}")
				


	


	if not math.isnan(RadRatio):
		RadBool = True

	else:
		url = r"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,tic_id,default_flag,pl_orbper,pl_ratdor,pl_ratror,pl_trandep,pl_radj,st_rad,pl_orbsmax+from+ps+where+tic_id=%27TIC%20"+TIC+"%27+and+default_flag=1+&format=csv"
		print(url)
		params = pd.read_table(url,sep=',')

		oneparams = params.iloc[0]

		print(oneparams)

		RadRatio = oneparams['pl_ratror']

	#with open(f'./{shortname}/SetupFiles/{shortname}.priors.1', 'a', encoding='utf-8') as file:
	with open(f'./SetupFiles/{shortname}.priors.1', 'a', encoding='utf-8') as file:
	#with open(f'{shortname}.priors.1', 'a', encoding='utf-8') as file:
		file.writelines(f"period_0 {periodEBLM} \n")
		file.writelines(f"tc_0 {tc_tic} \n")
		file.writelines(f"t14_0 {tc_dur} \n")
		if RadBool:
			file.writelines(f"p_0 {RadRatio} \n")


	#with open(f'./{shortname}/SetupFiles/{shortname}.priors.1', 'r') as file:
	with open(f'./SetupFiles/{shortname}.priors.1', 'r') as file:
	#with open(f'{shortname}.priors.1', 'r') as file:

		data = file.read()
		print(data)

	#print(f"Created /{shortname}/SetupFiles/{shortname}.priors.1")
	print(f"Created ./SetupFiles/{shortname}.priors.1")
	#os.chdir(f'./{shortname}')



	
if __name__ == "__main__":
    starval =0
    shortname = "HD209458"
    ticID = 420814525
    
    maxi=0
    #os.chdir('..')

    command =f"pwd"
    print("starting loc")
    os.system(command)

    #MakeNewPriorFile(starval,maxi,shortname,ticID)
    #UpdatePriorFile(shortname,ticID)

    ftm = pd.read_csv(r'Alison_new_ranked_SNRv3_best_copy.csv')
    tic = ftm.TIC
    
    for starval, ticID in enumerate(tic[:2]):
        shortname = ftm['Common_Name'][starval]
        shortname = shortname.replace(" ", "")
        os.chdir(f'./{shortname}')	
        UpdatePriorFile(shortname,ticID)
        os.chdir('..')


    	




