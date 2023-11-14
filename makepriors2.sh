#!/usr/bin/env bash
#SBATCH --time=00:30:00 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=11
#SBATCH --mem=16g
#SBATCH --job-name=priors
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
source ~/source.sh
ml idl
cd /home/duck.18/idl/EXOFASTv2/Benchmark 

idl -e "mkticsed, '420814525',priorfile='HD209458b.priors',sed='HD209458b.sed'" 
mv HD209458b.priors HD209458b/SetupFiles/HD209458b.priors 
mv HD209458b.sed HD209458b/Data/HD209458b.sed 
wait
