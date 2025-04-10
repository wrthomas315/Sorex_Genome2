#!/bin/bash

#environment : mask_genomes

#Sorex araneus
#first model repetitive regions using RepeatModeler
mkdir ../data/1_mask/mSorAra2/
gunzip ../data/0_refs/mSorAra2.pri.cur.20220707.fasta.gz
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -qy ../data/0_refs/mSorAra2.pri.cur.20220707.fasta -c 2 -o ../data/1_mask/mSorAra/mSorAra2_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -u <slurmusername> -qy ../data/0_refs/mSorAra2.pri.cur.20220707.fasta -c 2 -o ../data/1_mask/mSorAra/mSorAra2_masked
gzip ../data/0_refs/mSorAra2.pri.cur.20220707.fasta.gz
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/mSorAra2/1_masking2/RepeatMasker/mSorAra2_RM/mSorAra2.masked.fa ../data/1_mask/mSorAra2/mSorAra2.masked.fa
faToTwoBit ../data/1_mask/mSorAra2/mSorAra2.masked.fa  ../data/1_mask/mSorAra2/mSorAra2.masked.2bit
faSize -detailed -tab ../data/1_mask/mSorAra2/mSorAra2.masked.fa > ../data/1_mask/mSorAra2/mSorAra2.masked.chromsize


#Also run for Mustela nigripes
mkdir ../data/1_mask/musNig/
gunzip ../data/0_refs/GCF_022355385.1_MUSNIG.SB6536_genomic.fna.gz
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -qy ../data/0_refs/GCF_022355385.1_MUSNIG.SB6536_genomic.fna -c 2 -o ../data/1_mask/musNig/musNig_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -u <slurmusername> -qy ../data/0_refs/GCF_022355385.1_MUSNIG.SB6536_genomic.fna -c 2 -o ../data/1_mask/musNig/musNig_masked
gzip ../data/0_refs/GCF_022355385.1_MUSNIG.SB6536_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/musNig/1_masking2/RepeatMasker/GCF_022355385_RM/GCF_022355385.masked.fa ../data/1_mask/musNig/musNig.masked.fa
faToTwoBit ../data/1_mask/musNig/musNig.masked.fa  ../data/1_mask/musNig/musNig.masked.2bit
faSize -detailed -tab ../data/1_mask/musNig/musNig.masked.fa > ../data/1_mask/musNig/musNig.masked.chromsize


#And species for DESCHRAMBLER
#eriEur
mkdir ../data/1_mask/eriEur/
gunzip ../data/0_refs/GCA_950295315.1_mEriEur2.1_genomic.fna.gz
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -qy ../data/0_refs/GCA_950295315.1_mEriEur2.1_genomic.fna -c 2 -o ../data/1_mask/eriEur/$result'_masked'
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -u <slurmusername> -qy ../data/0_refs/GCA_950295315.1_mEriEur2.1_genomic.fna -c 2 -o ../data/1_mask/eriEur/$result'_masked'
gzip ../data/0_refs/GCA_950295315.1_mEriEur2.1_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/eriEur/1_masking2/RepeatMasker/GCA_950295315_RM/GCA_950295315.masked.fa ../data/1_mask/eriEur/eriEur.masked.fa
faToTwoBit ../data/1_mask/eriEur/eriEur.masked.fa  ../data/1_mask/eriEur/eriEur.masked.2bit
faSize -detailed -tab ../data/1_mask/eriEur/eriEur.masked.fa > ../data/1_mask/eriEur/eriEur.masked.chromsize

#phyDis
mkdir ../data/1_mask/phyDis/
gunzip ../data/0_refs/GCA_004126475.3_mPhyDis1.pri.v3_genomic.fna.gz
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -qy ../data/0_refs/GCA_004126475.3_mPhyDis1.pri.v3_genomic.fna -c 2 -o ../data/1_mask/phyDis/phyDis_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -u <slurmusername> -qy ../data/0_refs/GCA_004126475.3_mPhyDis1.pri.v3_genomic.fna -c 2 -o ../data/1_mask/phyDis/phyDis_masked
gzip ../data/0_refs/GCA_004126475.3_mPhyDis1.pri.v3_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/phyDis/1_masking2/RepeatMasker/GCA_004126475_RM/GCA_004126475.masked.fa ../data/1_mask/phyDis/phyDis.masked.fa
faToTwoBit ../data/1_mask/phyDis/phyDis.masked.fa  ../data/1_mask/phyDis/phyDis.masked.2bit
faSize -detailed -tab ../data/1_mask/phyDis/phyDis.masked.fa > ../data/1_mask/phyDis/phyDis.masked.chromsize

#rhiFer
mkdir ../data/1_mask/rhiFer/
gunzip ../data/0_refs/GCA_004115265.3_mRhiFer1_v1.p_genomic.fna.gz 
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -u <slurmusername> -qy ../data/0_refs/GCA_004115265.3_mRhiFer1_v1.p_genomic.fna -c 2 -o ../data/1_mask/rhiFer/rhiFer_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -qy ../data/0_refs/GCA_004115265.3_mRhiFer1_v1.p_genomic.fna -c 2 -o ../data/1_mask/rhiFer/rhiFer_masked
gzip ../data/0_refs/GCA_004115265.3_mRhiFer1_v1.p_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/rhiFer/1_masking2/RepeatMasker/GCA_004115265_RM/GCA_004115265.masked.fa ../data/1_mask/rhiFer/rhiFer.masked.fa
faToTwoBit ../data/1_mask/rhiFer/rhiFer.masked.fa  ../data/1_mask/rhiFer/rhiFer.masked.2bit
faSize -detailed -tab ../data/1_mask/rhiFer/rhiFer.masked.fa > ../data/1_mask/rhiFer/rhiFer.masked.chromsize

#talOcc
mkdir ../data/1_mask/talOcc/
gunzip ../data/0_refs/GCA_014898055.3_MPIMG_talOcc4v2.1_genomic.fna.gz
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -qy ../data/0_refs/GCA_014898055.3_MPIMG_talOcc4v2.1_genomic.fna -c 2 -o ../data/1_mask/talOcc/talOcc_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -u <slurmusername> -qy ../data/0_refs/GCA_014898055.3_MPIMG_talOcc4v2.1_genomic.fna -c 2 -o ../data/1_mask/talOcc/talOcc_masked
gzip ../data/0_refs/GCA_014898055.3_MPIMG_talOcc4v2.1_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/talOcc/1_masking2/RepeatMasker/GCA_014898055_RM/GCA_014898055.masked.fa ../data/1_mask/talOcc/talOcc.masked.fa
faToTwoBit ../data/1_mask/talOcc/talOcc.masked.fa  ../data/1_mask/talOcc/talOcc.masked.2bit
faSize -detailed -tab ../data/1_mask/talOcc/talOcc.masked.fa > ../data/1_mask/talOcc/talOcc.masked.chromsize

#sunEtr
mkdir ../data/1_mask/sunEtr/
gunzip ../data/0_refs/GCA_024139225.1_mSunEtr1.pri.cur_genomic.fna.gz
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -u <slurmusername> -qy ../data/0_refs/GCA_024139225.1_mSunEtr1.pri.cur_genomic.fna -c 2 -o ../data/1_mask/sunEtr/sunEtr_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -qy ../data/0_refs/GCA_024139225.1_mSunEtr1.pri.cur_genomic.fna -c 2 -o ../data/1_mask/sunEtr/sunEtr_masked
gzip ../data/0_refs/GCA_024139225.1_mSunEtr1.pri.cur_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/sunEtr/1_masking2/RepeatMasker/GCA_024139225_RM/GCA_024139225.masked.fa ../data/1_mask/sunEtr/sunEtr.masked.fa
faToTwoBit ../data/1_mask/sunEtr/sunEtr.masked.fa  ../data/1_mask/sunEtr/sunEtr.masked.2bit
faSize -detailed -tab ../data/1_mask/sunEtr/sunEtr.masked.fa > ../data/1_mask/sunEtr/sunEtr.masked.chromsize

#conCri
mkdir ../data/1_mask/conCri/
gunzip ../data/0_refs/GCA_000260355.1_ConCri1.0_genomic.fna.gz
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -qy ../data/0_refs/GCA_000260355.1_ConCri1.0_genomic.fna -c 2 -o ../data/1_mask/conCri/conCri_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -u <slurmusername> -qy ../data/0_refs/GCA_000260355.1_ConCri1.0_genomic.fna -c 2 -o ../data/1_mask/conCri/conCri_masked
gzip ../data/0_refs/GCA_000260355.1_ConCri1.0_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/conCri/1_masking2/RepeatMasker/GCA_000260355_RM/GCA_000260355.masked.fa ../data/1_mask/conCri/conCri.masked.fa
faToTwoBit ../data/1_mask/conCri/conCri.masked.fa  ../data/1_mask/conCri/conCri.masked.2bit
faSize -detailed -tab ../data/1_mask/conCri/conCri.masked.fa > ../data/1_mask/conCri/conCri.masked.chromsize

#galPyr
mkdir ../data/1_mask/galPyr/
gunzip ../data/0_refs/GCA_019455555.1_Gpyr_1.0_genomic.fna.gz   
python ./scripts_bin/JUSTmodel_genomes_toga_NEWUNITY.py -qy ../data/0_refs/GCA_019455555.1_Gpyr_1.0_genomic.fna -c 2 -o ../data/1_mask/galPyr/galPyr_masked
#then mask the genome with Repeat Masker
python ./scripts_bin/JUSTmask_genomes_toga.py -u <slurmusername> -qy ../data/0_refs/GCA_019455555.1_Gpyr_1.0_genomic.fna -c 2 -o ../data/1_mask/galPyr/galPyr_masked
gzip ../data/0_refs/GCA_019455555.1_Gpyr_1.0_genomic.fna
#and then turn these files to 2bit and get chromosome files for alignments
cp ../data/1_mask/galPyr/1_masking2/RepeatMasker/GCA_019455555_RM/GCA_019455555.masked.fa ../data/1_mask/galPyr/galPyr.masked.fa
faToTwoBit ../data/1_mask/galPyr/galPyr.masked.fa  ../data/1_mask/galPyr/galPyr.masked.2bit
faSize -detailed -tab ../data/1_mask/galPyr/galPyr.masked.fa > ../data/1_mask/galPyr/galPyr.masked.chromsize
