#!/bin/bash

#
#SBATCH --job-name=Files_WRT
#SBATCH --output=hmmclean.txt
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -p long-40core

#likely need to insall below or update paths for UNITY

module load shared
module load HmmCleaner/0.180750
module load anaconda/3
module load gnu-parallel/6.0

#uodate paths, also generic parralelize (not spelled right) parralelization of HmmCleaner

parallel --verbose -j 40 --link "HmmCleaner.pl {1}" ::: ./PATH/To/ALIGNMENTS/*.fasta
