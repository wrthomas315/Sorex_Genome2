#! /bin/bash

while read l; do
	echo "hyphy absrel --alignment ~/Sorex_Genome2/data/8_taper_align/alignments4hp_taper2/noFS"$l".macse.aligned.fa.taper.fa --tree ~/Sorex_Genome2/data/9_trees/unrootedtrees/noFS"$l".macse.aligned.fa.taper.fa.nh --branches All --output ~/Sorex_Genome2/analysis/3_hyphy/2_abs_dropout/jsons/"$l".ABSREL.explore.json" >>  ~/Sorex_Genome2/analysis/3_hyphy/2_abs_dropout/jobList 
done < ~/Sorex_Genome2/data/000_miscFiles/ENSTlist4hypy.txt

#then parallelize jobList however you would like across your cluster

while read l; do
	echo "hyphy absrel --alignment ~/Sorex_Genome2/data/8_taper_align/alignments4hp_taper2/noFS"$l".macse.aligned.fa.taper.fa --tree ~/Sorex_Genome2/data/9_trees/unrootedtrees/noFS"$l".macse.aligned.fa.taper.fa.nh --branches FORE --output ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/jsons/"$l".ABSREL.explore.json" >>  ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/jobList
done < ~/Sorex_Genome2/data/000_miscFiles/ENSTlist4hypy.txt
