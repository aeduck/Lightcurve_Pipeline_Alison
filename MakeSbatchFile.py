
import numpy as np
import pandas as pd
import sys
import os
import math

def MakeSbatch(shortname,pathexo = "/home/duck.18/idl/EXOFASTv2/Benchmark",email=""):
	"""
		Creates a Sbatch file for running EXOFASTv2.

		Args:
			shortname: The short name of the system.
			pathexo: string of the path where you want all your files to live (string)
			email: email adress for Unity to send emails to


		Returns:
			None.

		Produces:
			A starting sbatch file for less than 5 hours
	"""

	shortname_list = shortname.lower().rstrip()
	lower_shortname = shortname_list.replace("-", "")

	idlcommand = f'{lower_shortname}_1'

	boilerplate =["#!/usr/bin/env bash\n","\n","#SBATCH --time=5:00:00 \n","#SBATCH --nodes=1 \n","#SBATCH --ntasks-per-node=10\n","#SBATCH --mem=16g\n",f"#SBATCH --job-name=run{lower_shortname}_1\n","#SBATCH --mail-type=ALL\n",f"#SBATCH --mail-user={email}\n","\n","ml idl\n","source ~/source.sh\n"]

	
	with open(f'./SetupFiles/run{lower_shortname}_1.sh', 'w', encoding='utf-8') as file:

		for b in boilerplate:
			file.writelines(b)

		file.writelines(f"cd {pathexo}/{shortname}/SetupFiles\n")
		file.writelines(f'idl -e "{idlcommand}" \n')
		file.writelines('wait \n')
	print()


if __name__ == "__main__":
    shortname = "HD209458"
    ticID = 420814525

    MakeSbatch(shortname)