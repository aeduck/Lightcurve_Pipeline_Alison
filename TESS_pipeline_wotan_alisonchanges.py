#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sys
import astropy
import astropy.units as u
import lightkurve as lk
from lightkurve import search_targetpixelfile
import scipy as scipy
from scipy import stats
import matplotlib.pyplot as plt
import os
from astropy.time import Time
from io import StringIO
# only need these if wanting to look up values through code rather than entering them manually
import requests
from bs4 import BeautifulSoup as bsoup
import re # for parsing through ExoFOP (unconfirmed planets)
import wotan


# Search Lightkurve for your target first, then creat lightcurve files for the particular sector. Example at the end.
# general functions that can be used outside of making lcs

def bindata(x,y,binwidth):  
	
	'''
	This is a custom binning function. 
	x: time or phase array
	y: flux array
	binwidth: width of bins - I have NOT included astropy units here so make sure it's in the correct units!!!
	'''

	xbin = np.arange(np.nanmin(x),np.nanmax(x),binwidth)
	ybin = np.ones(len(xbin))
	#errbin = np.ones(len(xbin))
	
	for i in range(len(xbin)):
		mask = np.abs(x-xbin[i]) < binwidth/2
		#print(mask)
		ybin[i] = np.nanmedian(y[mask])

		#errbin[i] = np.sqrt(np.sum(err[mask]**2))
		
	# in case there are empty bins, exlcude these    
	nanmask = ~np.isnan(ybin)
	ybin = ybin[nanmask]
	xbin = xbin[nanmask]
	
	#print('nanmask: ',nanmask)
	return xbin,ybin#,errbin 

def calc_t14(Rp,Rstar,b,period,sma,inc):
	
	'''
	This returns the transit duration in hours, for use if it's not listed on NEA/EXOFOP etc.
	Astropy units ARE used here so make sure to include units on all inputs!

	
	Rp: planet radius
	Rstar: star radius
	b: impact parameter   --- if this isn't listed, set b = 0 to calculate longest possible transit duration
	Period should be self explanatory
	sma: semi-major axis
	inc: inclination (make sure to include deg or rad astropy units for input)
	'''
	
	C = (1.+(Rp/Rstar))**2 - b**2
	t14 = period/np.pi * np.arcsin(Rstar*C/sma/(np.sin(inc.to('rad'))).value)/u.rad
	
	return t14.to('hour')

def dilute(contam):
	
	'''
	Returns dilution term for a given TESS contamination ratio (listed at top left on EXOFOP).
	Include this in the prior file as: dilute_1 0.0 dilute_term
	'''
	
	return (contam/(1+contam))*0.1


# main function for getting lightcurve

def create_light_curve(target, author, sector, period=None, duration=None, tc=None, targetname=None, 
	exposure = None,multisector = False, save = False, plot=True, binlc=False, binwidth = None, auto=False,
	system = None, planet = None, qualityFlag = False, depth = None, published = True):
	
	'''
	Create light curve files - currently for TESS only. If only a single sector is needed, use 
	create_light_curve directly. If lightcurves for all sectors/cadences are needed, use multi_sector.
	
	target: name of target system, e.g. 'TIC ###'
	author: pipeline name, e.g. 'SPOC'  
	sector: TESS sector
	period: planet period in days. include astropy units when manually inputting.
	duration: transit duration in hours. include astropy units when manually inputting.
	tc: time of transit midpoint as BJD. don't include units.
	targetname: optional name used for naming files. if not specified default to target
	exposure time: observing cadence. do not need to specify, used for multisector or if you want specific LC
	multisector: flag for creating lcs for multiple sectors. 
	save: flag to save lightcurves to text files    
	plot: flag to show plots in notebook
	binlc: flag to bin the lightcurve(s)
	binwidth: width of bins used if binning
	auto: automatically search the NEA for values for period, tc and transit duration
	system: system name in K2-## format. only needed if using automatic search for values, has to match NEA.
	planet: planet letter. only needed if using automatic NEA search for values.
	quality flag: if true only uses Lightkurve data that has zero bad quality flags. good for removing major spikes
	depth:
	published: Boolean for whether or not the planet has been published. Reverts to ExoFOP if false.
	'''
	
	############## 
	# added - erica
	
	# If running for multiple sectors, specify exposure time in case of multiple cadences in one sector
	# Can also specify exposure time for a single lc

	#Finding and downloading target information and light curve
	targetinfo = lk.search_lightcurve(target, author = author, sector = sector, exptime=exposure)   
	   
	if exposure != None:
		targetinfo = targetinfo[np.where(targetinfo.table['t_exptime']==exposure)].download()
	
	# if no exposure time is specified download first lk entry by default
	if exposure == None:
		exposure =targetinfo.table['t_exptime'][0]
		if len(targetinfo) > 1:
			print('Warning! There are multiple entries for this sector. Only using first entry.')            
		targetinfo = targetinfo.download()
		
	targetinfo = targetinfo.remove_nans()
	
	if qualityFlag == True:
		# flag to only use data that has no bad data quality flags
		targetinfo = targetinfo[np.where(targetinfo.quality==0)]  
	
	# just checking if it's doing what i want
	#print('Sector: %s. Exposure time: %s seconds.' %(sector, int(exposure)))
	
	if author == 'QLP'or 'TESS-SPOC':
	   
		bjd = targetinfo.meta['BJDREFI']+targetinfo.meta['TSTART']
		date = Time(bjd,format='jd').fits[:10].replace('-','')
	
	else:
		date = targetinfo.meta['DATE'].replace('-','')
	
	# if target name isn't specified use input target by default
	if targetname == None:
		targetname = target
	
	if auto == True:
		period,tc, duration,depth = get_params(targetname,planet,published=published,target=target)
		
	
	#################
	
	#Determining Period, Transit Duration, and Tc
	period = period.to('day').value
	transit_duration = duration.to('day').value
	tc = tc
	print('period: ',period,'duration: ',duration,'tc: ',tc) 
	#Creating a buffer for transit plotting
	buffer = 0.010 * transit_duration
	#print('buffer: ',buffer)
	
	#Creating target time and flux information
	time = np.array(targetinfo.time.value + targetinfo.meta['BJDREFI'])
	flux = np.array(targetinfo.flux.value)/np.median(np.array(targetinfo.flux.value))
	
	#Phase folding the light curve to identify transit properties
	phase = (time-tc)/period-np.floor((time-tc)/period)
	gt5 = phase > 0.5
	phase[gt5] = phase[gt5]-1.0
	
	#Determining the size of the phase transit
	transit_size_phase = ((transit_duration)/period+(buffer)/period)*3

	#Indicating where transits are located within the light curve
	in_transit = np.where((phase>-transit_size_phase/2) & (phase<transit_size_phase/2))
	out_of_transit = np.where(~((phase>-transit_size_phase/2) & (phase<transit_size_phase/2)))


	phase_s = (time-(tc+period/2))/period-np.floor((time-(tc+period/2))/period)
	gt5_s = phase_s > 0.5
	phase_s[gt5_s] = phase_s[gt5_s]-1.0

	transit_s_size_phase = ((transit_duration)/period) - .2*(((transit_duration)/period))
	in_eclipse = np.where((phase_s>-transit_s_size_phase/2) & (phase_s<transit_s_size_phase/2))
	#print("length of in_eclipse array")
	#print(len(in_eclipse[0]))

	out_of_transit = np.where( ~((phase>-transit_size_phase/2) & (phase<transit_size_phase/2)) & ~((phase_s>-transit_s_size_phase/2) & (phase_s<transit_s_size_phase/2)))

	
	#Creating masks to exclude transits from the flattening process
	input_mask = np.ones_like(phase, dtype=bool)
	input_mask[in_transit] = False 
	
	#Flattening the light curve for any stellar variability
	#Alison changed!!!

	mask = wotan.transit_mask(time=time, period=period, duration=transit_duration+0.05*transit_duration, T0=tc)
	flat_flux, trend_lc = wotan.flatten(time, flux, mask = mask, window_length=3 * transit_duration, break_tolerance=0.2,  return_trend=True,verbose=True)
	
	binflag = 'unbinned'
	
	# If binning data, overwrite time, phase, flux and flat flux with binned versions
	if binwidth != None:
		print('Binning to :', binwidth)
		binflag = 'binned-%s' %str(binwidth).replace(' ','')
		binwidth = binwidth.to('day').value
		bintime, flux = bindata(time,flux,binwidth)
		bintime, flat_flux = bindata(time,flat_flux,binwidth)
		time = bintime
		phase = ((time)-tc)/period-np.floor(((time)-tc)/period)
		gt5 = phase > 0.5
		phase[gt5] = phase[gt5]-1.0
		
		in_transit = np.where((phase>-transit_size_phase/2) & (phase<transit_size_phase/2))
		out_of_transit = np.where(~((phase>-transit_size_phase/2) & (phase<transit_size_phase/2)))   
	
	#Separating the individual transits
	transit_times = []
	lc_transits = []
	array_names = []
	in_transit_n = []
	floored_time = np.floor(((time)-tc)/period)
	unique_epochs = np.unique(floored_time)
	
	#print('unique_epochs: ',unique_epochs)
	
	for val in unique_epochs:            #Attaching the actual transit times to their location indicators
		t_time = tc + (val*period)
		transit_times.append(t_time)
	#print('transit_times: ', transit_times)
	for j in transit_times:      #Only including transits within the observation timeframe
		if j > time[0] and j < time[-1]:
			lc_transits.append(j)   
	#print('lc_transits: ', lc_transits)
	for i in lc_transits:                  
		index = lc_transits.index(i)         #Creating an array for naming each available transit
		transit_name = 'in_transit_%s' % index
		array_names.append(transit_name)
	
		
		#Creating an array for where each transit can be found
		#in_transit_n.append(np.where((time>i-3/2*transit_duration-buffer) & (time<i+3/2*transit_duration+buffer)))  
		low_cut = (1.5*transit_duration-buffer)/24.
		hi_cut = (1.5*transit_duration+buffer)/24.
		
		in_transit_n.append(np.where((time>i-low_cut) & (time<i+hi_cut)))
	
	#print('array_names: ',array_names)
	#print('in_transit_n: ', in_transit_n)
	
	#Combining both the name array and location array into a comprehensive dictionary 
	# (a dataset for each transit)
	final_transits = dict(zip(array_names, in_transit_n))  
		
	# per-point error calculation - specific to TESS



	error_value_raw = scipy.stats.median_abs_deviation(flux[out_of_transit],nan_policy='omit')   #error per point
	print("not flat err:",error_value_raw)
	error_value_raw_std = scipy.stats.tstd(flux[out_of_transit])   #error per point
	#print("not flat err std:",error_value_raw)

	#error_value = scipy.stats.median_abs_deviation(flat_flux[out_of_transit],nan_policy='omit')/0.68   #error per point #Joey's group had this line!!!
	
	error_value = scipy.stats.median_abs_deviation(flat_flux[out_of_transit],nan_policy='omit')   #error per point
	print('per-point error value out of transit', error_value)
	error_value_std = scipy.stats.tstd(flat_flux[out_of_transit])   #error per point
	#print('per-point error value out of transit std', error_value)
	
	errors = np.array(np.full((len(time), 1), error_value, dtype=float))

	error_value_s = scipy.stats.median_abs_deviation(flat_flux[in_eclipse],nan_policy='omit')  #error per point
	#print('per-point error value in eclipse:', error_value_s)
	error_value_s_std = scipy.stats.tstd(flat_flux[in_eclipse])  #error per point
	#print('per-point error value in eclipse std:', error_value_s)
	
	# comparison of rms to transit depth
	if auto == True:
		rms = np.sqrt(np.sum((1.-flat_flux[out_of_transit])**2.)/len(flat_flux[out_of_transit]))
		print('Transith depth: %f, RMS: %.4f, Depth/RMS: %.3f' %(depth, rms,depth/rms))

	if plot == True:

		#Plotting the flattened, masked flux over the original flux to see the effects of cleaning the data
		fig1, ax1 = plt.subplots(figsize=(16, 8))
		plt.scatter(time, flux,label='Unflat Flux')
		plt.scatter(time, flat_flux, facecolors = 'tan', s = 100, alpha = 0.5, edgecolors='#000000',label='Flat Flux')
		plt.title("%s: Sector %s, %ss %s - Flattened, Masked Flux Vs. Original Flux" % \
				  (targetname,sector,int(exposure),binflag)) 
		plt.xlabel("Time [BTJD days]")
		plt.ylabel("Normalized Flux")
		plt.legend()

		#Plotting the light curve of only the transits
		fig2, ax2 = plt.subplots(figsize=(16, 8))
		plt.scatter(time[in_transit], flat_flux[in_transit])
		plt.title("%s: Sector %s, %ss %s - Individual Transit Flux" % (targetname,sector,int(exposure),binflag))
		plt.xlabel("Time - 2457000 [BTJD days]")
		plt.ylabel("Normalized Flux")
		#plt.ylim(0.992,1.01)
		#plt.xlim(2.4584e6+40.25,2.4584e6+41)

		#Plotting Individual Transits
		for key, value in final_transits.items():   
			newvalue = np.array(value)
			fig3, ax3 = plt.subplots(figsize=(16, 8))
			ax3.ticklabel_format(useOffset=False)
			plt.scatter((phase[newvalue[0,:]]*period), flat_flux[newvalue[0,:]])
			plt.title("%s: Sector %s, %ss %s - %s" % (targetname,sector,int(exposure),binflag, key))
			plt.xlabel("Time Since Transit Center (Days)")
			plt.ylabel("Normalized Flux")
			#print('transit times:',time[newvalue[0]])
			plt.show()
			plt.close()

		#Plotting the phase folded light curve to easily observe transit
		fig4, ax4 = plt.subplots(figsize=(16, 8))
		plt.scatter(phase[in_transit], flat_flux[in_transit])
		plt.title("%s: Sector %s, %ss %s - Phase Folded Light Curve" % (targetname,sector,int(exposure),binflag)) 
		plt.xlabel("Phase")
		plt.ylabel("Normalized Fluxf")
		#plt.xlim(min(phase[in_transit]),max(phase[in_transit]))
		#plt.ylim(0.99,1.01)

		fig5, ax5 = plt.subplots(figsize=(16, 8))
		plt.scatter(time, flux,label='Unflat Flux')
		plt.scatter(time, flat_flux, facecolors = 'tan', s = 100, alpha = 0.5, edgecolors='#000000',label='Flat Flux')
		plt.title("%s: Sector %s, %ss %s - Flattened, Masked Flux Vs. Original Flux" % \
				  (targetname,sector,int(exposure),binflag)) 
		plt.xlabel("Time [BTJD days]")
		plt.ylabel("Normalized Flux")
		plt.legend()
		
		

	# Saving data files of masked light curve, phase folded light curve, and light curve 
	# of the transits by themselves
	
	if save == True:
		
		########## CHANGE SAVE PATH HERE ##########
		
		
		startdir = os.getcwd()
		print("Saving Figures!")
		if os.path.exists('./Data'):
			os.chdir('./Data')
		else:
			os.mkdir('./Data')
			os.chdir('./Data')

		if plot == False:
			fig4, ax4 = plt.subplots(figsize=(16, 8))
			plt.scatter(phase, flat_flux,color='blue')
			plt.scatter(phase[in_transit], flat_flux[in_transit],color='blue')
			plt.xlim(-1*(transit_size_phase*2),(transit_size_phase*2))
			plt.title("%s: Sector %s, %ss %s - Phase Folded Light Curve" % (targetname,sector,int(exposure),binflag)) 
			plt.xlabel("Phase")
			plt.ylabel("Normalized Fluxf")

			fig5, ax5 = plt.subplots(figsize=(16, 8))
			plt.scatter(time, flux,label='Unflat Flux')
			plt.scatter(time, flat_flux, facecolors = 'tan', s = 100, alpha = 0.5, edgecolors='#000000',label='Flat Flux')
			plt.title("%s: Sector %s, %ss %s - Flattened, Masked Flux Vs. Original Flux" % \
					  (targetname,sector,int(exposure),binflag)) 
			plt.xlabel("Time [BTJD days]")
			plt.ylabel("Normalized Flux")
			plt.legend()

			fig6, ax6 = plt.subplots(figsize=(16, 8))
			plt.scatter(phase_s, flat_flux)
			plt.scatter(phase_s[in_eclipse], flat_flux[in_eclipse],color='red')
			plt.xlim(-1*(transit_size_phase*2),(transit_size_phase*2))
			plt.title("%s: Sector %s, %ss %s - Secondary Phase Folded Light Curve" % (targetname,sector,int(exposure),binflag)) 
			plt.xlabel("Phase")
			plt.ylabel("Normalized Fluxf")


			if os.path.exists('./Figs'):
				os.chdir('./Figs')
			else:
				os.mkdir('./Figs')
				os.chdir('./Figs')
			fig4.savefig('%s_S%s_phase_folded_%ss%s.png' % (targetname,sector,int(exposure),binflag),dpi=400, bbox_inches="tight",format='png',facecolor='white')
			fig5.savefig('%s_S%s_%ss_comparison.png' % (targetname,sector,int(exposure)), bbox_inches="tight",facecolor='white')
			fig6.savefig('%s_S%s_phase_folded_secondary_%ss%s.png' % (targetname,sector,int(exposure),binflag),dpi=400, bbox_inches="tight",format='png',facecolor='white')

			file1 = open(f'{targetname}_S{sector}_SummaryStats_{exposure}s.txt',"w")
			summary_stats= np.array([f'Sector: {sector} Exposure time: {exposure} seconds \n', f'period: {period :.3f} duration: {duration :.3f} tc: {tc :.4f} \n',f'not flat per-point error value: {error_value_raw}\n',f'not flat per-point error value (std): {error_value_raw_std}\n',f'per-point error value: {error_value}\n',f'per-point error value (std): {error_value_std}\n',f'per-point error value in eclipse: {error_value_s} \n',f'per-point error value (std) in eclipse: {error_value_s_std} \n'])
			file1.writelines(summary_stats)
			file1.close()

			os.chdir('..')


		

		#fig4.savefig('%s_S%s_phase_folded_%ss%s.png' % (targetname,sector,int(exposure),binflag),dpi=400, bbox_inches="tight",format='png',facecolor='white')
		#fig5.savefig('%s_S%s_%ss_comparison.png' % (targetname,sector,int(exposure)), bbox_inches="tight",facecolor='white')

		#fig4.close()
		#ig5.close()
		plt.clf()
		plt.close('all')


		np.savetxt('n%s.TESS.TESS.%s.S%s.%ss%s.FullNotFlat.dat' % (date,targetname, sector,int(exposure),binflag), \
				   np.c_[time, flux, errors], delimiter=' ') 
		np.savetxt('n%s.TESS.TESS.%s.S%s.%ss%s.FullFlat.dat' % (date,targetname, sector,int(exposure),binflag), \
				   np.c_[time, flat_flux, errors], delimiter=' ') 
		np.savetxt('n%s.TESS.TESS.%s.S%s.%ss%s.SlimFlat.dat' % (date,targetname, sector,int(exposure),binflag), \
				   np.c_[time[in_transit], flat_flux[in_transit], errors[in_transit]], delimiter=' ')   
		
		# save individual transits to separate files

		'''
		for key, value in final_transits.items():
			newvalue = np.array(value)
			np.savetxt('n%s.TESS.TESS.%sSlimFlat.S%s.%ss%s.%s.dat' % (date,targetname, sector,int(exposure),binflag, key), \
					   np.c_[(phase[newvalue[0,:]]*period), flat_flux[newvalue[0,:]]], delimiter=' ') 
		'''
		os.chdir(startdir)





def GrabSector(name,ticID,sec,exp):
	
	

	ticName = "TIC "+str(ticID)


	lc = lk.search_lightcurve(ticName,author='SPOC', mission="TESS", sector=sec, exptime=exp)#J0525-55
	#print(lc)




	short_mission = [m[-2:] for m in lc.mission]
	


	url = 'https://exofop.ipac.caltech.edu/tess/download_planet.php?id=' + str(ticID)
	print(url)


	url = f'https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/{ticID}/info/?tce=1&sector={sec}'
	print(url)
	params = pd.read_json(url,orient='index')

	duration_1sig = (params['TDUR'][1] + params['INDUR'][1])/24 #duration +ingress
	tc_tic = params['TEPOCH'][1]
	tic_period = params['TPERIOD'][1]

	tc_tic  = float(params["TEPOCH"][1])+float(params["BJDREFI"][1]) # in JD
	#print(tc_tic)



	#print("Sector "+m)
	create_light_curve(ticName,'SPOC',sec,auto=False, exposure=exp,  period=tic_period*u.day, duration= duration_1sig*u.day, tc=tc_tic,
						save = True,plot=False,qualityFlag = True,system=name,planet=' ')
	print()



def missions(t):

	targetinfo = lk.search_lightcurve(f'TIC {t}', author='SPOC', mission="TESS")
	if not targetinfo:
		n = [np.nan, np.nan]
		return n
	print(targetinfo)
	mis = pd.DataFrame(np.array(targetinfo.table['mission','exptime']))
	n = mis.groupby('mission').min()


	sectors=[]
	for thing in n.index:
		sectors.append(thing[12:])
	

	exptimes=[]
	for thing in n.values:
		exptimes.append(thing[0])
	
	num_of_sectors = len(n)


	#counts = n.value_counts()
	#exptimes = [t[0] for t in counts.index]
	#fullinfo = np.transpose([exptimes,counts])
	fullinfo = np.transpose([sectors,exptimes])

	#print("fullinfo!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	#print(fullinfo)

	

	return [num_of_sectors,fullinfo]



def startTess(ticID,name,all120s=False):

	print()
	print()
	print()

	print("name")
	print(name)


	n = missions(ticID)
	print("Number of Sectors:", n[0])
	#print("n")
	#print(n)

	for ns in n[1]:
		sec = ns[0]
		exp = int(float(ns[1]))
		if int(sec) >= 69:
			print("As of October 2023 these systems have not been processed yet.")
			continue

		if all120s:
			exp = 120
			print(f"Sector: {sec} \nExptime: {120}")
			GrabSector(name,ticID,sec,exp)

		if not all120s:
			print(f"Sector: {sec} \nExptime: {exp}")
			GrabSector(name,ticID,sec,exp)



if __name__ == "__main__":
	name = "HD209458b"
	ticID = 420814525
	
	#name="WASP33b"
	#ticID=129979528
	#name="WASP-100b"
	#ticID=38846515
	tic=[ticID]


	#name="WASP-187b"
	#name="WASP-187b"
	#ticID=15692883
	#tic=[ticID]


	#GrabSector(name,ticID)
	#idl -e "mkticsed, '420814525',priorfile='HD209458.priors',sed='HD209458.sed'" 


	#ftm = pd.read_csv(r'FinalTargetsWithProperMotions.csv')
	#tic = ftm.TIC


	for t in tic:

		n = missions(t)
		print("Number of Sectors:", n[0])
		#print("n")
		#print(n)

		for ns in n[1]:
			sec = ns[0]
			exp = int(float(ns[1]))
			print(f"Sector: {sec} \nExptime: {exp}")
			GrabSector(name,ticID,sec,exp)