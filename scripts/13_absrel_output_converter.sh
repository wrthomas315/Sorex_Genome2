#!/bin/bash


#Foreground
ls ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/jsons/ > ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/jsonList.txt

while read j; do
	ENST=$(echo $j | awk -F'.' '{print $1}')
	sorUNCORR=$(awk '/sorex/{flag=1; next} flag && /Uncorrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j) 
	sorCORR=$(awk '/sorex/{flag=1; next} flag && /Corrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
	sunUNCORR=$(awk '/suncus/{flag=1; next} flag && /Uncorrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
	sunCORR=$(awk '/suncus/{flag=1; next} flag && /Corrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
	furUNCORR=$(awk '/putorius/{flag=1; next} flag && /Uncorrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
	furCORR=$(awk '/putorius/{flag=1; next} flag && /Corrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
	ermUNCORR=$(awk '/erminea/{flag=1; next} flag && /Uncorrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
	ermCORR=$(awk '/erminea/{flag=1; next} flag && /Corrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
	echo $ENST "sorex_araneus" $sorUNCORR $sorCORR >> ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/output.txt
	echo $ENST "suncus_etruscus" $sunUNCORR $sunCORR >> ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/output.txt
	echo $ENST "mustela_erminea" $ermUNCORR $ermCORR >> ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/output.txt
	echo $ENST "mustela_putorius_furo" $furUNCORR $furCORR >> ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/output.txt
done < ~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/jsonList.txt


#Drop out
ls ~/Sorex_Genome2/analysis/3_hyphy/2_abs_dropout/jsons/ > ~/Sorex_Genome2/analysis/3_hyphy/2_abs_dropout/jsonList.txt

while read j; do
	ENST=$(echo $j | awk -F'.' '{print $1}')
	cat jsons/$j | grep "test" | grep -v "results" | grep -v "tested" | grep -v "branch" | awk -F':' '{print $1}' | awk -F'"' '{print $2}' > specs.txt
	while read m; do
		sorUNCORR=$(awk -v key=$m '$0 ~ key {flag=1; next} flag && /Uncorrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j) 
		sorCORR=$(awk -v key=$m '$0 ~ key {flag=1; next} flag && /Corrected/{match($2, /[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?/); printf "%.10f\n", substr($2, RSTART, RLENGTH); flag=0; exit}' jsons/$j)
		echo $ENST $m $sorUNCORR $sorCORR >> ~/Sorex_Genome2/analysis/3_hyphy/2_abs_dropout/output.txt
	done < specs.txt
done < ~/Sorex_Genome2/analysis/3_hyphy/2_abs_dropout/jsonList.txt
