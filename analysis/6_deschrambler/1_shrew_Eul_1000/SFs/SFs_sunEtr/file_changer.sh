#!/bin/bash

sed 's/SUPER_//g' Orthology.Blocks > Orthology.Blocks.mod
sed 's/_//g' Orthology.Blocks.mod > Orthology.Blocks.mod0
COUNTER=0
while read j; do
	WOUNTER=$(( COUNTER ))
	COUNTER=$(( COUNTER + 1 ))
	echo $WOUNTER
	echo $COUNTER
	old=$(echo $j | awk '{print $1}' | awk -F"_" '{print $1$2}')
	echo $old
	new=$(echo $j | awk '{print $2}')
	sed "s/$old/$new/g" Orthology.Blocks.mod"$WOUNTER" > Orthology.Blocks.mod"$COUNTER"
done <sun.chrom.key

#cat APCF_sorAra.merged.map | grep . | grep "-" | paste -d ' ' - - | awk -F. -v OFS=; '{print $1 OFS $2}'
#for APCFmap
#cat APCF_sorAra.merged.map | grep . | grep "-" | paste -d ' ' - - | awk -F. -v OFS=";" '{print $2 OFS $3 OFS $4}' | awk -F: -v OFS=";" '{print $1 OFS $2 OFS $3}' | awk -F" " -v OFS=";" '{print $1 OFS $3 OFS $4}' | awk -F"SUPER_" -v OFS=";" '{print $1$2}' | sed 's/\+\;/1/g' | sed 's/\-\;/PLACE/g' | awk -F"-" -v OFS=";" '{print $1 OFS $2 OFS $3}' | sed 's/X/12/g' | sed 's/PLACE/-1/g' | awk -F";" -v OFS='\t' '{print $1 OFS $2 OFS $3 OFS $5 OFS $6 OFS $7 OFS $8 OFS $4}' > formatted_syntenyBlocks.txt
#for species to species comparison
#cat Orthology.Blocks_con | grep . | grep "-" | paste -d ' ' - - | awk -F. -v OFS=";" '{print $2 OFS $3 OFS $4}' | awk -F: -v OFS=";" '{print $1 OFS $2 OFS $3}' | awk -F" " -v OFS=";" '{print $1 OFS $3 OFS $4}' | awk -F"SUPER_" -v OFS=";" '{print $1$2}' | sed 's/\+\;/1/g' | sed 's/\-\;/PLACE/g' | awk -F"-" -v OFS=";" '{print $1 OFS $2 OFS $3}' | sed 's/PLACE/-1/g' | awk -F";" -v OFS='\t' '{print $1 OFS $2 OFS $3 OFS $5 OFS $6 OFS $7 OFS $8 OFS $4}' > sorcon_syntenyBlocks.txt
