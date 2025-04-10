#!/bin/bash

#run with slurm
#make sure Toga is installed https://github.com/hillerlab/TOGA
#need to update four nextflow configs, can change --cb to maximize your specific clusters memory

#run Toga
~/bin/TOGA/toga.py ~/data/4_clean_genome_align/mSorAra2/hg38.mSorAra2.masked.all.rr.clean.chain ~/data/000_miscFiles/hg38.v35.for_toga.bed ~/data/0_refs/hg38.2bit ~/data/1_mask/mSorAra2/mSorAra2.masked.2bit \
        --kt --pd ~/data/5_TOGA/mSorAra2/toga --nc ~/bin/TOGA/nextflow_config_files/  --cb 10,100,268 -i ~/data/000_miscFiles/hg38.v35.for_toga.isoforms.tsv \
         --u12 ~/data/000_miscFiles/hg38.U12sites.tsv --ms --cjn 21
#run for musNig as well
