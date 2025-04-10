##########################################################INFO###############################################################################################
#  mask_genomes.py			
#
#  Setup script for mask genomes with Repeat modeler and Repeat masker.
#
#  Author: Diana Moreno edited for use on Unity in Feb 2022. Bill Thomas edited 2022 for Unity paths, then edited for local use with indirect paths.
#  v.1.0 November 2019
#  v.2.0 January 2020 --> add xsmall and lib as Flags.
#  v.3.0 April 2021 --> TOGA
#  v.4.0 WRT indirect paths
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
#           -u <slurm username>
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
REPMODOUT = (MAINDIR + '/1_masking2/RepeatModeler/RM_output/consensi.fa.classified') #Assign the output from Repeat Modeler to a path for further analysis.

#3)Define input arguments:

def get_args():
	parser = argparse.ArgumentParser(description='Input needed for variant call annotation', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-qy', '--query', type=str, help='Path/to/query/assembly.fa with variant call', required=True)
	parser.add_argument('-u', '--user', type=str, help='Slurm username', required=True)
	parser.add_argument('-c', '--cores', type=int, help='Define number of cores to use for multi-threaded process', default=1, required=True) 
	parser.add_argument('-o', '--outdir', type=str, help='Define the path/to/output/directory', default='./') 
	args = parser.parse_args()
	QRY = args.query 
	CORES = args.cores
	OUTDIR = args.outdir
	USER = args.user
	return QRY, CORES, OUTDIR, USER

QRY, CORES, OUTDIR, USER = get_args() 

QUERY = os.path.abspath(QRY) #Set the complete path for query fasta file. Due to during all the scripts we are switching directories, by assigning the full path we avoid errors. 
PREFIX = os.path.basename(QRY).split(".")[0]

#4)Sanity check

print("The query genome is {}\n".format(QRY))
print("The query genome is located in {}\n".format(QUERY))
print("Output files will be generated with {} prefix".format(PREFIX))
print("For this run you are going to use {} cores".format(CORES))
print("Results will be saved on {}".format(OUTDIR))


#II. REPEAT MASKER -----------------------------------------------------------------

#  **IIa.- Run RepeatMasker using the denovo library built with Repeat Modeler.
#        - IMPORTANT: This script uses a python program. You need to make sure python is loaded and ready to go by activating conda.
#        -The following script was developed by David A. Ray

os.chdir(MAINDIR) #Moving from /1_Repeat_modeler directory to /2_Repeat_masker.
os.chdir('1_masking2')
os.mkdir('RepeatMasker')
os.chdir('RepeatMasker')

RM2BED=('rm2bed.sh')

with open(RM2BED, 'w') as f:
    f.write('#!/bin/bash' + '\n')
    f.write('#SBATCH --job-name=' + PREFIX + '_RM2' + '\n')
    f.write('#SBATCH --output=%x.%j.out' + '\n')
    f.write('#SBATCH --error=%x.%j.err' + '\n')
    f.write('#SBATCH --partition=cpu-long' + '\n')
    f.write('#SBATCH --cpus-per-task=16' + '\n')
    f.write('#SBATCH --nodes=1' + '\n') 
    f.write('#SBATCH --ntasks=1' + '\n')
    f.write('\n') 
    f.write('echo makes variables used in this script' + '\n')
    f.write('GENOME=' + PREFIX + '\n')
    f.write('USER=' + USER + '\n')
    f.write('RUNTYPE=' + PREFIX + '_RM' + '\n')
    f.write('DIR=' + MAINDIR + '/1_masking2/RepeatMasker/$RUNTYPE' + '\n')
    f.write('cd $DIR' + '\n')
    f.write('#Runs a python script RM2bed to generate one complete .bed file and several subfiles subdivided by TE class. Merges overlapping hits based using lower_divergence criterion.' + '\n')
    f.write('[ ! -f ${GENOME}_rm.bed ] && python RM2bed.py -d . -sp class -p ${GENOME} -o lower_divergence ${GENOME}.fa.align.gz' + '\n')
    f.write('module load gcc/9.3.0 bedtools/2' + '\n')
    f.write('RMPATH=../../bin/envs/mask_genomes/bin/' + '\n')
    f.write('gunzip -c ${GENOME}.fa.out.gz > ${GENOME}.fa.out' + '\n')
    f.write('perl $RMPATH/rmOutToGFF3.pl ${GENOME}.fa.out > ${GENOME}.gff' + '\n')
    f.write('bedtools maskfasta -soft -fi ' + QRY + ' -bed ${GENOME}.gff -fo ${GENOME}.masked.fa' + '\n')
    f.write( 'echo "Masking genome done!"' + '\n')

RMSLURM=('Repeatmasker_slurm.sh')



with open(RMSLURM, 'w') as f:
    f.write('#!/bin/bash' + '\n')
    f.write('#SBATCH --job-name=' + PREFIX + '_repmask' + '\n')
    f.write('#SBATCH --output=%x.%j.out' + '\n')
    f.write('#SBATCH --error=%x.%j.err' + '\n')
    f.write('#SBATCH --partition=cpu-long' + '\n')
    f.write('#SBATCH --nodes=1' + '\n')
    f.write('#SBATCH --ntasks=1' + '\n')
    f.write('#SBATCH --cpus-per-task=16' + '\n')
    f.write('\n')
    f.write('echo makes variables used in this script' + '\n')
    f.write('GENOME=' + QRY + '\n' + '\n')
    f.write('USER=' + USER + '\n')
    f.write('RUNTYPE=' + PREFIX + '_RM' + '\n')
    f.write('DIR=' + MAINDIR + '/1_masking2/RepeatMasker/$RUNTYPE' + '\n')
    f.write('mkdir -p $DIR' + '\n')
    f.write('cd $DIR' + '\n' + '\n')
    f.write('ln -s ' + QUERY + '\n' + '\n')
    f.write('echo "fa to 2bit"'+ '\n')
    f.write(FATOBIT + ' ' + QRY + ' ' + PREFIX + '.2bit' + '\n' + '\n')
    f.write('echo "Use slurm_clusterrunC8.py to generate all batches needed to run RepeatMasker and the doLift.sh to complile the results."'+ '\n')
    f.write('python3 ' + CLUSTERRUN + ' -i ' + QRY + ' -b 50 -lib ' + REPMODOUT + ' -dir . -xsmall' + '\n')
    f.write('echo "Submits the RepeatMasker jobs created by slurm_clusterrunC8.py"'+ '\n')
    f.write('bash qsub.sh'+ '\n')
    f.write('#Creates a list of jobIDs to keep track of the RepeatMasker batches being run.' + '\n')
    f.write('jobIDs=""; for i in `squeue  | grep ' + PREFIX + """ | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs=\"${jobIDs:1}\" """ + '\n')
    f.write('#Submits the doLift script to the queue but holds it until all jobs in the jobIDs list (the RepeatMasker batches) have cleared.' + '\n')
    f.write('sleep 10m'+ '\n')
    f.write('X=`squeue -u $USER | wc -l`  '+ '\n')
    f.write('while [ $X -gt 2 ]; do sleep 1m; X=`squeue -u $USER | wc -l`; done'+ '\n')
    f.write('sbatch doLift.sh' + '\n')
    f.write('#Creates a list of job IDs to keep track of the doLift script.' + '\n')
    f.write('jobIDs=""; for i in `squeue  | grep ' + PREFIX + """ | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs=\"${jobIDs:1}\" """ + '\n')
    f.write('#Submits the rm2bed script to the queue but holds it until all jobs in the jobIDs list (the doLift job) have cleared.' + '\n')
    f.write('sleep 10m'+ '\n')
    f.write('X=`squeue -u $USER | wc -l`  '+ '\n')
    f.write('while [ $X -gt 2 ]; do sleep 1m; X=`squeue -u $USER | wc -l`; done'+ '\n')
    f.write('sbatch ' + MAINDIR + '/1_masking2/RepeatMasker/rm2bed.sh' + '\n')
  
   
subprocess.check_call (['sbatch', 'Repeatmasker_slurm.sh']) 

