#!/bin/bash

#Get files for each found in () from https://genome.senckenberg.de//download/TOGA/
cd ~/Sorex_Genome2/data/7_alignments
#get transcript list for each species and make input directory file
rm ./input_dirs
while read k; do
	cd ~/Sorex_Genome2/data/7_alignments/
	cd $k
	awk '{print $2}' orthology_classification.tsv | sort -u > genelist
	pwd >> ../input_dirs
done < CG_speclistPLUS
#grab only unique transcripts
cd ~/Sorex_Genome2/data/7_alignments
cat ./*/genelist | sort -u > ./completeList
cp ./completeList ./completeList.tmp
sed '$ d' ./completeList.tmp > ./completeList
rm -f ./completeList.tmp


#Create alignment jobList
#need input dirs - list of species input directories
#completeList -
#might need to update line MACSE_TWO = "java -jar ~/Sorex_Genome2/bin/macse_v2.06.jar"... in ~/Sorex_Genome2/scripts/scripts_bin/extract_codon_alignmentMACMOD.py if macse in different location
while read l; do
        echo "python ~/Sorex_Genome2/scripts/scripts_bin/extract_codon_alignmentMACMOD.py input_dirs ~/Sorex_Genome2/data/000_miscFiles/hg38.v35.for_toga.bed "$l" --reference_2bit ~/Sorex_Genome2/data/0_refs/hg38.2bit -d --max_copies 1 -o ~/Sorex_Genome2/data/7_alignments/"$l".macse.aligned.fa --temp_dir ~/Sorex_Genome2/data/7_alignments/tmp" >>  jobList done < ./completeList


#then parallelize on your own depending on your cluster settings

#once run, select longest transcripts for each gene to run hyphy on
#make ENSG list you will loop through
awk '{print $1}' ~/Sorex_Genome2/data/000_miscFiles/hg38.v35.for_toga.isoforms.tsv | tail -n +2 | sort -u  > ~/Sorex_Genome2/data/7_alignments/ENSG_list
rm LongestTranscript_full
rm LongestTranscript_4PAML
#
while read j; do
	rm ENSTs
        #get a list of ENSTs associated with each ENSG
        cat ~/Sorex_Genome2/data/000_miscFiles/hg38.v35.for_toga.isoforms.tsv | grep  $j | awk '{print $2}' > ENSTs
        #see how  long this key is  and save  as  variable
	ENSTfilesize=$(cat  ENSTs | wc -l)
	#if 0
	if [[ "$ENSTfilesize" -lt 1 ]]; then
		echo -e $j"\tno_alignment"  >> LongestTranscript_full
	#if 1
	elif [[ "$ENSTfilesize" -lt 2 ]]; then
		f=$(cat ENSTs)
		speciesCount=$(cat "~/Sorex_Genome2/data/7_alignments/"$f".macse.aligned.fa" | wc -l)
		if [[ "$speciesCount" -gt 58 ]]; then
			echo -e $j"\t"$f >> LongestTranscript_full
			echo -e $j"\t"$f >> LongestTranscript_4PAML
		elif [[ "$speciesCount" -lt 59 ]]; then
			echo -e $j"\tnot_75per"  >> LongestTranscript_full
		fi
	#if greater than 1
	elif [[ "$ENSTfilesize" -gt 1 ]]; then
		#set a counter
		counter=0
		rm longest_intermediate
		#create new while looking at ENSTs rather than ENSG
		while read k; do
			p=$(cat "~/Sorex_Genome2/data/7_alignments/"$k".macse.aligned.fa" | grep "REF" -A1 | grep -o "[a-zA-Z]*"| wc -c)
			if [[ "$p" -gt "$counter"  ]]; then
				counter=$p
                        	#if higher than counter, overwrite final j file with new longest
                        	echo $k > longest_intermediate
                	fi
		done < ENSTs
		l=$(cat longest_intermediate)
		m=$(cat "~/Sorex_Genome2/data/7_alignments/"$l".macse.aligned.fa" | wc -l)
		if [[ "$m" -gt 58 ]]; then
			echo -e $j"\t"$l >> LongestTranscript_full
			echo -e $j"\t"$l >> LongestTranscript_4PAML
		elif [[ "$m" -lt 59 ]]; then
			echo -e $j"\tnot_75perL"  >> LongestTranscript_full
		fi
	fi
done < ENSG_list
