#!/bin/bash

#Download sequencing from SRA into ~/Sorex_Genome2/data/00_transcriptomics/seq/

#Then use id files to trim adapters and remove low quality reads
#!/bin/bash

#Get only ids
awk '{print $1}' ~/Sorex_Genome2/data/00_transcriptomics/ids/ids_cortex.txt | grep "HFS" > ~/Sorex_Genome2/data/00_transcriptomics/ids/ids_fastKal_c.txt
awk '{print $1}' ~/Sorex_Genome2/data/00_transcriptomics/ids/ids_hippocampus.txt | grep "HFS" > ~/Sorex_Genome2/data/00_transcriptomics/ids/ids_fastKal_h.txt

#make directory for trimmed
mkdir ~/Sorex_Genome2/data/00_transcriptomics/seq/trimmed
#loop through both with fastp (cortex shown)
while read j; do
        echo $j
	fastp --in1 "~/Sorex_Genome2/data/00_transcriptomics/seq/"$j"_R1_001.fastq" --in2 "~/Sorex_Genome2/data/00_transcriptomics/seq/"$j"_R2_001.fastq" --out1 "~/Sorex_Genome2/data/00_transcriptomics/seq/trimmed/"$j"_1.trimmed.fastq" --out2 "~/Sorex_Genome2/data/00_transcriptomics/seq/trimmed/"$j"_2.trimmed.fastq" --cut_tail --html $j".html" --json $j".json" 2> $j".log"
done <~/Sorex_Genome2/data/00_transcriptomics/ids/ids_fastKal_c.txt

#Make transcriptome index
REF=~/Sorex_Genome2/data/0_refs/GCF_027595985.1_mSorAra2.pri_rna.fna
IDX=~/Sorex_Genome2/data/0_refs/GCF_027595985.1_mSorAra2.pri_rna.idx
kallisto index -i $IDX $REF

#Pseudo align to transcriptome using kallisto
while read j; do
        echo $j
        kallisto quant -i "~/Sorex_Genome2/data/0_refs/GCF_027595985.1_mSorAra2.pri_rna.idx" -o "~/Sorex_Genome2/analysis/5_expression/Cortex/TranscriptAbundances/"$j "~/Sorex_Genome2/data/00_transcriptomics/seq/trimmed/"$j"_1.trimmed.fastq" "~/Sorex_Genome2/data/00_transcriptomics/seq/trimmed/"$j"_2.trimmed.fastq" > "~/Sorex_Genome2/analysis/5_expression/Cortex/TranscriptAbundances/"$j".kalistolog.txt"
done <~/Sorex_Genome2/data/00_transcriptomics/ids/ids_fastKal_c.txt
