#!/bin/bash

#Get genes needed to run meme on
awk '{print $2}' shrewtotal.txt | grep -v "X" > memeList

#make directory for trees need
mkdir ~/Sorex_Genome2/data/9_trees/unrootedtrees4meme/

#label trees with species for each gene
while read j; do
	cat ~/Sorex_Genome2/analysis/3_hyphy/hyphytotal.txt | grep $j | awk '{print $3}' > specs
	cp ~/Sorex_Genome2/data/9_trees/unrootedtrees4hp/noFS$j".macse.aligned.fa.taper.fa.nh" ~/Sorex_Genome2/data/9_trees/unrootedtrees4meme/
	while read m; do
		sed -i "s/$m/$m{FORE}/g" ~/Sorex_Genome2/data/9_trees/unrootedtrees4meme/noFS$j".macse.aligned.fa.taper.fa.nh"
	done < specs
done < ./memeList

#create meme job list
while read l; do
	echo "hyphy meme --alignment ~/Sorex_Genome2/data/8_taper_align/alignments4hp_taper2/noFS"$l".macse.aligned.fa.taper.fa --tree ~/Sorex_Genome2/data/9_trees/unrootedtrees4meme/noFS"$l".macse.aligned.fa.taper.fa.nh --branches FORE --output ~/Sorex_Genome2/analysis/4_meme/jsons/"$l".ABSREL.meme.json" >>  jobList
done < ./memeList 

#parallelize and run, then move on to sorting jsons for information
