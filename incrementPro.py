#Read in pro file and increment
#example: python incrementpro KELT15 1 n
#That would create fitKELT15_2.pro
import os
import sys
import re
import fileinput


def checkInput(message, error, goalString=False,goalInt=False,goalFloat=False):
	"""
		Gets input from the user and checks that it is valid.

		Args:
			message: The message to be displayed to the user.
			error: The error message to be displayed if the input is invalid.
			goalString: Whether the input should be a string.
			goalInt: Whether the input should be an integer.
			goalFloat: Whether the input should be a float.

		Returns:
			The valid input from the user.
	"""
	
	badInput = True
	counter = 0
	while badInput:
		value = input(message)
		if goalString:
			if any(c.isalpha() for c in value):
				badInput = False
				return value
	
		elif goalInt:
			
			if value.isdigit():
				if float(value)%1 == 0:
					value = int(value)
					badInput = False
					return value
					
		elif goalFloat:
			try:
				value = float(value)
			except:
				badInput=True
			else:
				if isinstance(value,(float)):
					badInput = False
					return value

		if badInput:
			print(error)
			print()
			counter = counter +1

		if counter > 5:
			print("I give up!")
			quit()
		

def increment(system=False,currentVal=False,sbatch=False,sbatch_time=False):
	"""
		Increments the run number and creates a new PRO file and Sbatch file.

		Args:
			system: The name of the system. Must be  a string.
			currentVal: The current run number. Must be an integer
			sbatch: Whether or not to create an Sbatch file using 'y' or 'n'. Must be a string.
			sbatch_time: The time limit for the Sbatch file in hours. Must be a positive number.

		Returns:
			None.

		Produces:
			system_nextnumber.pro, system.prior.nextnumber, and runsystem_nextnumber.sh
  	"""

	if isinstance(system,bool):
		system = checkInput("Enter the name of your system: ","Error: That was not a system name",goalString = True)

	if isinstance(currentVal,bool):
		currentVal = checkInput("Enter the integer of your completed run: ","Error: That was not an integer",goalInt = True)
	
	if isinstance(sbatch,bool):
		sbatch = checkInput("Would you like an sbatch file? (y/n): ","Error: That was not y or n",goalString = True)

	if sbatch == 'y' or 'yes' in sbatch.lower():
		sbatch_q = True
	elif sbatch == 'n' or 'no' in sbatch.lower():
		sbatch_q = False

	if sbatch_q:
		if isinstance(sbatch_time,bool):
			sbatch_time = checkInput("How long do you want your run in hours: ","Error: That was not a number",goalFloat = True)


	nextVal = currentVal + 1

	shortname_list = system.lower().rstrip()
	lower_shortname = shortname_list.replace("-", "")

	print(f"Making {lower_shortname}_{nextVal}.pro, creating run{lower_shortname}_{nextVal}.sh,  and copying {system}.priors.{nextVal}")

	#Copying the priors.NextNumber file from the fitresults to my set up files folder

	command =f"cp ./{system}/Results/fit{system}_{currentVal}/{system}.priors.{nextVal} ./{system}/SetupFiles/{system}.priors.{nextVal}"
	os.system(command)


	#Makeing the next pro file in my SetUpFiles folder
	

	with open(f'./{system}/SetupFiles/{lower_shortname}_{currentVal}.pro', 'r') as file:
		data = file.read().replace(f'{lower_shortname}_{currentVal}', f'{lower_shortname}_{nextVal}')
		data = data.replace(f'{system}.priors.{currentVal}',f'{system}.priors.{nextVal}')
		data = data.replace(f"path+'../Results/fit{system}_{currentVal}' + path_sep()+'fit{system}_{currentVal}.',",f"path+'../Results/fit{system}_{nextVal}' + path_sep()+'fit{system}_{nextVal}.',")
		if sbatch_q:
			sbatch_time_seconds = int(sbatch_time*60*60)
			data = re.sub('maxtime = .*?,',f'maxtime = {sbatch_time_seconds},',data, flags=re.DOTALL) #actually keeps rest of line!
			data = re.sub('maxtime=.*?,',f'maxtime = {sbatch_time_seconds},',data, flags=re.DOTALL) #actually keeps rest of line!



	with open(f'./{system}/SetupFiles/{lower_shortname}_{nextVal}.pro', 'w') as file:
		file.write(data)


	if sbatch_q:
		with open(f'./{system}/SetupFiles/run{lower_shortname}_{currentVal}.sh', 'r') as file:
			data2 = file.read()
			data2 = re.sub(f"#SBATCH --time=.*?:00:00",f"#SBATCH --time={str(int(sbatch_time+4))}:00:00",data2) # this format will get rid of rest of line
			data2 = data2.replace(f'{lower_shortname}_{currentVal}', f'{lower_shortname}_{nextVal}')

			
		with open(f'./{system}/SetupFiles/run{lower_shortname}_{nextVal}.sh', 'w') as file:
			file.write(data2)



if __name__ == "__main__":	
	increment()