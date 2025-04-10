#!/bin/bash
#SBATCH --output=batch_%A_%a.out
#SBATCH --error=batch_%A.%a.err
#SBATCH --job-name=doBlast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu-long

#Note: all in same script here but would run separately
#Note2: Worked on Unity, future user will have to change sbatch configs and configs (login1) below to own cluster settings
#Note3: If process breaks down and have failed jobs, rerun failed jobs in psl folder and use -continue flag from where it broke down (ex. -continue chainMerge)
#Note4: Make sure ./axtChain/run/chain output doesn't have empty files, would suggest a "part" failed along the way

#For TOGA S. araneus and M. nivalis
mkdir ../data/3_genome_align/sorAra
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/sorAra/
cd ../data/3_genome_align/sorAra/ 
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_sor -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

mkdir ../data/3_genome_align/musNig
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/musNig/
cd ../data/3_genome_align/musNig/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_mus -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

#For aligning Eulipotyphla genomes to Sorex for chromosome evolution
mkdir ../data/3_genome_align/eriEur
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/eriEur/
cd ../data/3_genome_align/eriEur/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_eri -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

mkdir ../data/3_genome_align/phyDis
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/phyDis/
cd ../data/3_genome_align/phyDis/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_phy -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

mkdir ../data/3_genome_align/rhiFer
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/rhiFer/
cd ../data/3_genome_align/rhiFer/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_rhi -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

mkdir ../data/3_genome_align/talOcc
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/talOcc/
cd ../data/3_genome_align/talOcc/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_tal -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

mkdir ../data/3_genome_align/sunEtr
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/sunEtr/
cd ../data/3_genome_align/sunEtr/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_sun -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

mkdir ../data/3_genome_align/conCri
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/conCri/
cd ../data/3_genome_align/conCri/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_con -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log

mkdir ../data/3_genome_align/galPyr
cp ./scripts_bin/doBlastzChain_WRT_git.pl ../data/3_genome_align/galPyr/
cd ../data/3_genome_align/galPyr/
./doBlastzChain_WRT_git.pl ../2_DEF/DEF_gal -verbose=2 -noDbNameCheck -workhorse=login1 -bigClusterHub=login1 -skipDownload -dbHost=login1 -smallClusterHub=login1 -trackHub -fileServer=login1 -syntenicNet > do.log
