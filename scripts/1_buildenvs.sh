#!/bin/bash

#load your clusters version of anaconda
module load anaconda

#get ucsc kent resources
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ../bin/

#build a masking environment and a toga environment
conda env create -f ../bin/envs/mask_genomes.yml
conda env create -f ../bin/envs/toga.yml

