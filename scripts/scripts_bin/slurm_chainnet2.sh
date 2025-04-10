#!/bin/bash
#SBATCH --job-name=chain_cleaner
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=196:00:00
#SBATCH -p cpu-long


#export PATH="/home/tlama_umass_edu/bin:$PATH"
export PATH="/work/tlama_umass_edu/bin:$PATH"

REFERENCE=/home/tlama_umass_edu/project_bat1k_longevity/data/hg38/hg38.fa
GENOME=/home/tlama_umass_edu/project_bat1k_longevity/data/myoNig/mNig.masked.fa
REF2=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/hg38.2bit
GEN2=/project/tlama_umass_edu/projects/project_shrew_genome/data/genomes/hg38/trackData/myoVel/myoVel_masked.2bit
RFOUT=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/trackData/myoVel/cleaning/hg38.myoVel_masked.all.rr.chain
REFX=$(basename ${REFERENCE} | awk -F '.' '{print $1}')
QUEX=$(basename ${GENOME} | awk -F '.' '{print $1}')
Rsize=/work/tlama_umass_edu/project_shrew_genome/data/genomes/hg38/hg38.chrom.sizes
Qsize=/project/tlama_umass_edu/projects/project_shrew_genome/data/genomes/hg38/trackData/myoVel/myoVel.masked.chromsize
#Get genome size

#faSize -detailed ${GENOME} > ${QUEX}.size
#faSize -detailed ${REFERENCE} > ${REFX}.size

#Run chainCleanner
#chainCleaner ${RFOUT} ${REF2} ${GEN2} ${QUEX}_clean.chain ${QUEX}_clean.bed -tSizes=${Rsize} -qSizes=${Qsize} -linearGap=medium

#Run ChainNet
#chainNet output_clean.chain ${REFX}.size ${QUEX}.size ${REFX}.net ${QUEX}.net -rescore -tNibDir=${REF2} -qNibDir=${GEN2} -linearGap=medium
#chainNet ${RFOUT} ${Rsize} ${Qsize} ${REFX}.net ${QUEX}.net

#chainCleaner ${RFOUT} ${REF2} ${GEN2} ${QUEX}_clean.chain ${QUEX}_clean.bed -net=hg38.net  linearGap=medium
#chainCleaner ${RFOUT} ${REF2} ${GEN2} ${QUEX}_clean.chain ${QUEX}_clean.bed -tSizes=${Rsize} -qSizes=${Qsize} -linearGap=medium
chainNet ${RFOUT} ${Rsize} ${Qsize} ${REFX}.net ${QUEX}.net
#NetFilter

#NetFilterNonNested.perl -doUCSCSynFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 ref.query.net.gz > ref.query.filtered.net
