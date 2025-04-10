#!/bin/bash

#
#SBATCH --job-name=cleaner
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=336:00:00
#SBATCH -p cpu-long

REFERENCE=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/hg38.2bit
GENOME=/project/tlama_umass_edu/projects/project_shrew_genome/data/genomes/hg38/trackData/hypMon/hypMon_masked.2bit
CHAINOUT=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/trackData/hypMon/cleaning/hg38.hypMon_masked.all.chain
RFOUT=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/trackData/hypMon/cleaning/hg38.hypMon_masked.all.rr.chain

python /project/tlama_umass_edu/bin/GenomeTools/GenomeAlignmentTools/src/RepeatFiller.py -c ${CHAINOUT} -T2 ${REFERENCE} -Q2 ${GENOME} -o ${RFOUT}
