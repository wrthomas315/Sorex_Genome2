#!/bin/bash

#
#SBATCH --job-name=cleaner
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=336:00:00
#SBATCH -p cpu-long

#need lastz if did not align on own http://www.bx.psu.edu/~rsharris/lastz/
#get genome alignment tools, instructions here https://github.com/hillerlab/GenomeAlignmentTools and set up in ../bin

#make a directory for clean genome alignments for each species
mkdir ../data/4_clean_genome_align/sorAra

#set references and output, below for shrew but do for every species needed (Mustela nigripes and Eulipotyphla species)
REFERENCE=../data/0_refs/hg38.2bit
GENOME=../data/1_mask/mSorAra2/mSorAra2.masked.2bit
CHAINOUT=../data/3_genome_align/mSorAra2/hg38.mSorAra2.masked.all.chain
RFOUT=../data/4_clean_genome_align/mSorAra2/hg38.mSorAra2.masked.all.rr.chain
REFX=$(basename ${REFERENCE} | awk -F '.' '{print $1}')
QUEX=$(basename ${GENOME} | awk -F '.' '{print $1}')
Rsize=../data/0_refs/hg38.chrom.sizes
Qsize=../data/1_mask/mSorAra2/mSorAra2.masked.chromsize

ql
#Repeat fill
python ./bin/GenomeTools/GenomeAlignmentTools/src/RepeatFiller.py -c ${CHAINOUT} -T2 ${REFERENCE} -Q2 ${GENOME} -o ${RFOUT}

#Build nets
chainNet ${RFOUT} ${Rsize} ${Qsize} ${REFX}.net ${QUEX}.net

#Clean chains
chainCleaner ${RFOUT} ${REFERENCE} ${GENOME} ${QUEX}_clean.chain ${QUEX}_clean.bed -net=hg38.net  linearGap=medium

#for deschambler species, build nets again with cleaned chains
chainNet ${QUEX}_clean.chain ${Rsize} ${Qsize} ${REFX}.post.net ${QUEX}.post.net
