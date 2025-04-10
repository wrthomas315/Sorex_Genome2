#!/bin/bash

#
#SBATCH --job-name=chainer
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=336:00:00
#SBATCH -p cpu-long

AXTOUT=/home/tlama_umass_edu/GenomeTools/myoTest/time_reduced/myoNig_hg38TR_alignment2.axt
REFERENCE=/home/tlama_umass_edu/project_bat1k_longevity/data/annotations/hg38/hg38.fa
QUERY=/home/tlama_umass_edu/project_bat1k_longevity/data/annotations/myoNig/mNig.masked.fa
CHAINOUT=/home/tlama_umass_edu/GenomeTools/myoTest/time_reduced/myoNig_hg38TR_output.chain

axtChain -linearGap=medium ${AXTOUT} -faT ${REFERENCE} -faQ ${QUERY} ${CHAINOUT}
