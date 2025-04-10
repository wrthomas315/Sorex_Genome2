#!/bin/bash

#identify every json
ls ~/Sorex_Genome2/analysis/4_meme/jsons/ > jsonList.txt
rm MEMEoutput.txt

#loop through everyone
while read j; do
	#get gene associated an species
	ENST=$(echo $j | awk -F'.' '{print $1}')
	cat hyphytotal.txt | grep $ENST | awk '{print $3}' > specs.txt
	#loop through each and find the species and the EBF
	while read m; do
		rm nums.txt
		awk -v spec="$m" '{
			if(index($0, spec) && count == 0) {
				count = 1;
				next;
			}
			if(index($0, spec) && count == 1) {
				count = 2;
				next;
			}
			if(count == 1) {
				print;
			}
		}' ~/Sorex_Genome2/analysis/4_meme/jsons/$j | grep "EBF" | awk '{print $3 FS $5}' | awk '{if (match($0, /^([0-9]+).*:([0-9]+(\.[0-9]+)?(e\+[0-9]+)?),/, arr)) {gsub(/1e\+[0-9]+/, "10000000000000000000000000000", arr[2]); if (arr[2] > 100) print arr[1], arr[2]}}' > nums.txt
		#print out the Gene EBF species
		while read v; do
			echo $ENST $m $v >> MEMEoutput.txt
		done < nums.txt
	done < specs.txt
done < jsonList.txt

#Now have output for every species
#create a key that lumps genes and sites
awk '{key=$1 FS $3; if (!(key in seen)) seen[key]=$2; else seen[key]=seen[key]","$2; count[key]++} END {for (key in seen) if (count[key] >= 2) print key, seen[key]}' MEMEoutput.txt > convergencesiteTotal.txt
#now one that lumps genes and sites with just those including shrew
awk '{key=$1 FS $3; if (!(key in seen)) seen[key]=$2; else seen[key]=seen[key]","$2; count[key]++} END {for (key in seen) if (count[key] >= 2) print key, seen[key]}' MEMEoutput.txt |  grep "sorex" > convergencesiteShrew.txt#

#Now have sites for making figures
