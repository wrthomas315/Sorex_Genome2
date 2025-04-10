#!/bin/bash

#do for Sorex araneus for selection tests, and both Mustela nigripes and Sorex aranaues for BUSCO
#first make a document of all the genes said to be found  in species
awk '{print $4}'  ~/Sorex_Genome2/data/5_TOGA/mSoraAra2/query_gene_spans.bed > geneIDs
mkdir js
mkdir finalJs
#read this file as you  are going to loop through the genes
while read j; do
	#first create a variable that has a tab then E, so say reg_11 does not come up for a search of reg_1
        r=$j"	E"
	#now search for this gene in the orthology document and save results to a file
	cat ../orthology_classification.tsv | grep "$r" > ./js/"$j"
	#looper=$(cat ./js/"$r" | wc -l)
	#take only the 4th column, which has the ENST transcripts for each gene
	awk '{print $4}' ./js/"$j" > ./js/x_"$j"
	#start a counter, that as you loop through each transcript, the longest will be saved
	counter=0
	#so open this file and l oop through the ENST
	for i in $(cat ./js/x_"$j")
	do
		#search for each ENST in the query file, look at the next line, count alphabetic characters
		p=$(cat ~/Sorex_Genome2/data/5_TOGA/mSoraAra2/query_prot.fasta | grep "$i" -A1 | grep -o "[a-wy-zA-WY-Z]*"| wc -c)
		#first result will beat counter of 0, each after will have to beat new counter
		if [[ $p -gt $counter  ]]
		then
			counter=$p
			#if higher than counter, overwrite final j file with new longest
			cat ~/Sorex_Genome2/data/5_TOGA/mSoraAra2/query_prot.fasta | grep -w "$i" -A1 > ./finalJs/$j
		fi
	done
done < geneIDs

#combine the genes
while read j; do
        cat ./finalJs/$j  >>  LongestTranscript.fasta
done < geneIDs

#remove temp folders
rm -rd ./finalJs
rm -rd ./js
