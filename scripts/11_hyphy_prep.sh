#!/bin/bash
while read l; do
	for i in {1..40}
		do
		sed -n "$i"p ~/Sorex_Genome2/data/000_miscFiles/spec_key > spec
		re=$(awk '{print $1}' spec)
		#echo $re
		place=$(awk '{print $2}' spec)
		#echo $place
		sed -i "s/.*>$re.*/>$place/" ~/Sorex_Genome2/data/8_taper_align/noFS$l".macse.aligned.fa.taper.fa"
		done
done  < ~/Sorex_Genome2/data/000_miscFiles/ENSTlist4hypy.txt
