
##########################################################INFO###############################################################################################
#  mask_genomes.py			
#
#  Setup script for mask genomes with Repeat modeler and Repeat masker.
#
#  Author: Diana Moreno edited for use on Unity in Feb 2022, edited by Bill Thomas 2022
#  v.1.0 November 2019
#  v.2.0 January 2020 --> add xsmall and lib as Flags.
#  v.3.0 April 2021 --> TOGA
#  v.4.0 WRT for use on Unity and local for github
#
#  Description:
#
#	This script will perform the following steps for TOGA functional annotation:
#		I.- Identify de novo repeats with RepeatModeler.
#        The output for this step will be stored at 1_Repeat_modeler directory
#	   II.- Mask  genomes with RepeatMasker, using the RepeatModeler output.
#  
#    syntax: python variant_call_annotation.py
#           -qy <path/to/query/assembly.fa> \
#           -c <# number of cores to use> \
################################################################################################################################################################



#----------------------------------------------------------------- SETTING ENVIRONMENT -----------------------------------------------------------------

#1)Import modules:

import argparse
import os
import os.path
import subprocess
from subprocess import check_output
import sys 
import time


#2)Set path to software packages:
FATOBIT = ('../../bin/faToTwoBit')
CLUSTERRUN = ('./slurm_clusterrunC8.py')
REPMODE= ('../../bin/envs/mask_genomes/bin/') 
REPMASK = ('../../bin/envs/mask_genomes/bin/RepeatMasker') 
AUGUSTUS = ('../../bin/envs/mask_genomes/bin/augustus')
MAINDIR = os.getcwd() 


#3)Define input arguments:

def get_args():
	parser = argparse.ArgumentParser(description='Input needed for variant call annotation', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-qy', '--query', type=str, help='Path/to/query/assembly.fa with variant call', required=True)
	parser.add_argument('-c', '--cores', type=int, help='Define number of cores to use for multi-threaded process', default=1, required=True) 
	parser.add_argument('-o', '--outdir', type=str, help='Define the path/to/output/directory', default='./') 
	args = parser.parse_args()
	QRY = args.query 
	CORES = args.cores
	OUTDIR = args.outdir
	return QRY, CORES, OUTDIR

QRY, CORES, OUTDIR = get_args() 

QUERY = os.path.abspath(QRY) #Set the complete path for query fasta file. Due to during all the scripts we are switching directories, by assigning the full path we avoid errors. 
PREFIX = os.path.basename(QRY).split(".")[0]

#4)Sanity check

print("The query genome is {}\n".format(QRY))
print("The query genome is located in {}\n".format(QUERY))
print("Output files will be generated with {} prefix".format(PREFIX))
print("For this run you are going to use {} cores".format(CORES))
print("Results will be saved on {}".format(OUTDIR))

#5)Create working directories for each task
print ('Creating main working directory')
os.mkdir('1_masking2')
os.chdir('1_masking2')

#***I. RepeatModeler***

os.mkdir('RepeatModeler')
os.chdir('RepeatModeler')

#  **Ia.- Create qsub for Repeatmodeler

REPMODENAME=(PREFIX + '_RM_slurm.sh')

with open(REPMODENAME, 'w') as f:
        f.write('#!/bin/bash' + '\n')
        f.write('#SBATCH --job-name=' + PREFIX + '_RM' + '\n')
        f.write('#SBATCH --output=%x.%j.out' + '\n')
        f.write('#SBATCH --error=%x.%j.err' + '\n')
        f.write('#SBATCH --partition=cpu-long' + '\n')
        f.write('#SBATCH --nodes=2' + '\n')
        f.write('#SBATCH --cpus-per-task=16' + '\n')
        f.write('#SBATCH --ntasks=1' + '\n')
        f.write('#SBATCH --time=196:00:00' + '\n')
        f.write('\n' + '\n')
        f.write('module load gcc/9.3.0 ncbi-rmblastn/2.9.0ls' + '\n')
        f.write('echo "Build database using"' + QUERY + 'as a genome reference' + '\n')
        f.write((REPMODE) + '/BuildDatabase -name ' + PREFIX + '_repmod ' + '-engine ncbi ' + QUERY + '\n')
        f.write('echo "Running repeat modeler..."' + '\n' + '\n')
        f.write((REPMODE) + '/RepeatModeler -database ' + PREFIX + '_repmod -engine ncbi -pa ' + str(CORES) + '\n')
        f.write('mv RM* RM_output' + '\n')
        f.write('rm -r RM_output/round-*' + '\n')
        f.write('echo "Repeat modeler finished on $ (date)" > repmod_done.txt' + '\n')
        
#  **Ib.-  Excecute sbatch to run RepeatModeler and safety check points.

#print("Running repeat modeler, this might take a while...")

subprocess.check_call(['sbatch', REPMODENAME])

print("Running repeat modeler, this might take a while...")

while not os.path.exists('repmod_done.txt'):  #Whit this loop Python is paused until RepeatModeler finish.
        time.sleep(1)
        if os.path.exists('repmod_done.txt'):
                print ('Repeat Modeler finished')

#  **IIIb.-  Check if the output file required for RepeatMasker exists, if not the script will stop.
#    NOTE: This code is going to work only if you are located in the RepeatModeler directory.

try:
        open('./RM_output/consensi.fa.classified')
except IOError:  #file not found
        print ('File not found, Repeat Modeler failed. \n This Python script will die!')
        sys.exit()

REPMODOUT = (MAINDIR + '/1_masking2/RepeatModeler/RM_output/consensi.fa.classified') #Assign the output from Repeat Modeler to a path for further analysis.

