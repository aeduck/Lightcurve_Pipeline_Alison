import numpy as np
import pandas as pd
import sys
import os
from astroquery.eso import Eso
from astropy import table
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy import units as u
from pyvo.dal import tap
import requests
import cgi
import json
import getpass

def downloadURL(file_url, dirname='.', filename=None, session=None):
    """Method to download a file, either anonymously (no session or session not "tokenized"), or authenticated (if session with token is provided).
       It returns: http status, and filepath on disk (if successful)"""
    if not  os.path.exists('./archive'):
            os.mkdir('./archive')
    if dirname != None:
        if not os.access(dirname, os.W_OK):
            print("ERROR: Provided directory (%s) is not writable" % (dirname))
            sys.exit(1)
      
    if session!=None:
        response = session.get(file_url, stream=True)
    else:
        # no session -> no authentication
        response = requests.get(file_url, stream=True)

    # If not provided, define the filename from the response header
    if filename == None:
        contentdisposition = response.headers.get('Content-Disposition')
        if contentdisposition != None:
            value, params = cgi.parse_header(contentdisposition)
            filename = params["filename"]

        # if the response header does not provide a name, derive a name from the URL
        if filename == None:
            # last chance: get anything after the last '/'
            filename = file_url[file_url.rindex('/')+1:]

    # define the file path where the file is going to be stored
    if dirname == None:
        filepath = filename
    else:
        filepath = dirname + '/' + filename

    if response.status_code == 200:
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=50000):
                f.write(chunk)

    return (response.status_code, filepath)



def getDispositionFilename( response ):
	"""Get the filename from the Content-Disposition in the response's http header"""
	contentdisposition = response.headers.get('Content-Disposition')

	if contentdisposition == None:
		return None
	value, params = cgi.parse_header(contentdisposition)
	#print(params)
	filename = params["filename"]
	return filename

def writeFile( response ):
	"""Write on disk the retrieved file"""
	if response.status_code == 200:
		# The ESO filename can be found in the response header
		filename = getDispositionFilename( response )
		if not  os.path.exists('./archive'):
			os.mkdir('./archive')

		# Let's write on disk the downloaded FITS spectrum using the ESO filename:
		with open("archive/"+filename, 'wb') as f:
			f.write(response.content)
		return filename 

def ESOdownload(ticID,shortname,inst):


	ESO_TAP_OBS = "http://archive.eso.org/tap_obs"

	tapobs = tap.TAPService(ESO_TAP_OBS)

	#AlisonDuck
	#jsK2GsNKQWBU8h

	# Defining the position via SESAME name resolver, and the search radius
	#target = "EBLM J0525-55"

	target = shortname

	
	pos = SkyCoord.from_name(target)
	# pos now contains the coordinates of NGC 4666
	print("SESAME coordinates for %s: %s (truncated to millidegrees)\n" % (target, pos.to_string()))

	RA_val = pos.ra.degree 
	DEC_val = pos.dec.degree
	

	sr = 1/60. # search radius of 1 arcmin, always expressed in degrees

	# Cone search: looking for footprints of reduced datasets intersecting a circle of 1' around J0525-55

	query = f"""SELECT *
	FROM ivoa.ObsCore
	WHERE intersects(s_region, circle('', %f, %f, %f))=1
	AND obs_collection = '{inst}' """ % (RA_val , DEC_val, sr)

	print(query)



	res = tapobs.search(query=query, maxrec=1000)
	print("Num matching datasets: %d" % (len(res)))

	
	command =f"pwd"
	os.system(command)

	print(res['dp_id'])
	dp_ids = res['dp_id']
	command =f"pwd"
	os.system(command)

	for dp in dp_ids:

		status_code, filepath = downloadURL(f'http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{dp}','./archive')
		
		with open(filepath, 'r') as file:
			data = file.readlines()
			ancillary = False
			for i,line in enumerate(data):
				if 'ANCILLARY.HARP' in line:
					ancillary = True
				elif 'ANCILLARY.CCF' in line:
					ancillary = True
				if ancillary:
					if 'ivo://eso.org/ID?' in line:
						new_tar = line
						break
		os.system(f'rm {filepath}')


		#print("--------------------------------------------")
		#print(new_tar)
		#print(new_tar[75:-6])
		dp = new_tar[75:-6]
		file_url = f'https://dataportal.eso.org/dataportal_new/file/{dp}'
		

		response = requests.get(file_url)
		filename = writeFile( response )
		if filename:
			print("Saved file: %s" % (filename))
		else:
			print("Could not get file (status: %d)" % (response.status_code))




if __name__ == "__main__":

	shortname = "HD209458"
	ticID = 420814525
	
	inst='ESPRESSO'
	ESOdownload(ticID,shortname,inst)

	inst='HARPS'
	ESOdownload(ticID,shortname,inst)
