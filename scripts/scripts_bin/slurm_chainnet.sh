#!/bin/bash
#SBATCH --job-name=chain_cleaner
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=196:00:00
#SBATCH -p gpu-long

CHAINOUT=chain_output.chain
REFERENCE=/home/tlama_umass_edu/project_bat1k_longevity/data/hg38/hg38.fa
GENOME=/home/tlama_umass_edu/project_bat1k_longevity/data/myoNig/mNig.masked.fa
REF2=/home/tlama_umass_edu/project_bat1k_longevity/data/annotations/hg38/hg38.2bit
GEN2=/home/tlama_umass_edu/project_bat1k_longevity/data/annotations/myoNig/myoNig_masked.2bit
RFOUT=/home/tlama_umass_edu/GenomeTools/myoTest/time_reduced/myoNig_hg38TR_refill.chain
REFX=$(basename ${REFERENCE} | awk -F '.' '{print $1}')
QUEX=$(basename ${GENOME} | awk -F '.' '{print $1}')
Rsize=/home/tlama_umass_edu/project_bat1k_longevity/data/annotations/hg38/hg38.chrom.sizes
Qsize=/home/tlama_umass_edu/project_bat1k_longevity/data/annotations/myoNig/myoNig.size
#Get genome size

#faSize -detailed ${GENOME} > ${QUEX}.size
#faSize -detailed ${REFERENCE} > ${REFX}.size

#Run chainCleanner
#chainCleaner ${RFOUT} ${REF2} ${GEN2} ${QUEX}_clean.chain ${QUEX}_clean.bed -tSizes=${Rsize} -qSizes=${Qsize} -linearGap=medium

#Run ChainNet
#chainNet output_clean.chain ${REFX}.size ${QUEX}.size ${REFX}.net ${QUEX}.net -rescore -tNibDir=${REF2} -qNibDir=${GEN2} -linearGap=medium
chainNet ${RFOUT} ${Rsize} ${Qsize} ${REFX}.net ${QUEX}.net
#NetFilter

#NetFilterNonNested.perl -doUCSCSynFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 ref.query.net.gz > ref.query.filtered.net
