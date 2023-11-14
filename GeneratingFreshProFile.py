
import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
#boilerplate =["#!/usr/bin/env bash\n","\n","#SBATCH --time=5:00:00 \n","#SBATCH --nodes=1 \n","#SBATCH --ntasks-per-node=10\n","#SBATCH --mem=16g\n",f"#SBATCH --job-name=run{lower_shortname}_1\n","#SBATCH --mail-type=ALL\n","#SBATCH --mail-user=duck.18@osu.edu\n","\n","ml idl\n","source ~/source.sh\n"]


def MakeProFile(ticID,shortname,dirpath='',templatefile ='../templatepro.txt'):
    """
        Creates a PRO file for EXOFASTv2.

        Args:
            ticID: The TIC ID of the system.
            shortname: The name of the system.
            subdirpath: string subdirectories for exofast file path
                        example format f"['Benchmark','{shortname}','SetupFiles']"

        Returns:
            None.
        Produces:
            An EXOFASTv2 pro file for our system. By default it includes the SED and applys the MIST isochrones.
     """
    #template =["pro Q_Y, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin\n","if n_elements(maxsteps) eq 0 then maxsteps = 15000\n","if n_elements(nthin) eq 0 then nthin=50\n", f"path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['Benchmark','{shortname}','SetupFiles'])\n","\n","exofastv2, nplanets=1, tranpath=path+'../Data/*FullFlat.dat',","rvpath=path+'../Data/*.rv',","$\n","\tpriorfile=path+'Z.priors.Y', $\n","\tprefix=path+'../Results/fitZ_Y' + path_sep()+'fitZ_Y.', maxsteps=maxsteps, $\n","\tnthin=nthin, nthreads = 10, ntemps=8,maxtime=10800, circular=[0],fitthermal=['TESS'], fittran=[1]",", fitrv=[1]",", $\n",f"\tdebug=debug, verbose=verbose, mistsedfile = path+'../Data/{shortname}.sed\n","end\n"]
    if dirpath == '':
        dirpath = f"['Benchmark','{shortname}','SetupFiles']"
    
    

    if not templatefile:
        boilerplate =["pro Q_Y, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin\n","if n_elements(maxsteps) eq 0 then maxsteps = 15000\n","if n_elements(nthin) eq 0 then nthin=50\n", f"path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir={dirpath})\n","\n","exofastv2, nplanets=1, tranpath=path+'../Data/*FullFlat.dat',","rvpath=path+'../Data/*.rv',","$\n","\tpriorfile=path+'Z.priors.Y', $\n","\tprefix=path+'../Results/fitZ_Y' + path_sep()+'fitZ_Y.', maxsteps=maxsteps, $\n","\tnthin=nthin, nthreads = 10, ntemps=8,maxtime=10800, circular=[0],fitthermal=['TESS'], fitreflect=['TESS'], fittran=[1]",", fitrv=[1]",", $\n",f"\tdebug=debug, verbose=verbose, mistsedfile = path+'../Data/{shortname}.sed\'\n","end\n"]
    else:
        with open(templatefile, 'r') as tempfile:
            boilerplate = tempfile.readlines()
        print(boilerplate)
    


    #print(shortname)
    #print(os.system("pwd"))
    shortname_list = shortname.lower().rstrip()
    lower_shortname = shortname_list.replace("-", "")

    if not os.path.exists(f'./SetupFiles'):
        print("made folder!")
        os.mkdir('./SetupFiles')

    with open(f'./SetupFiles/{lower_shortname}_1.pro', 'w', encoding='utf-8') as file:
    #with open(f'{lower_shortname}_1.pro', 'w', encoding='utf-8') as file:
        for b in boilerplate:
            b = b.replace('{dirpath}', f'{lower_shortname}')
            b = b.replace('{shortname}', f'{shortname}')
            b = b.replace('Q', f'{lower_shortname}')
            b = b.replace('Z', f'{shortname}')
            b = b.replace('Y', '1')
            if 'rvpath' in b or 'fitrv' in b:
                if os.path.exists('./archive'):
                    file.writelines(b)
            else:
                file.writelines(b)
                
            

    print(f"Created SetUpFiles/{lower_shortname}_1.pro")
    print()
    #print(f"{lower_shortname}_1.pro")


if __name__ == "__main__":
    shortname = "HD209458b"
    ticID = 420814525
    MakeProFile(ticID,shortname)