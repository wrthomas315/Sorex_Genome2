#!/bin/bash

busco --in ~/Sorex_Genome2/data/6_longestTranscript/mSorAra2/LongestTranscript.fasta  --out ~/Sorex_Genome2/analysis/1_busco/mSorAra2/ -l mammalia_odb10 -m prot  -f
busco --in ~/Sorex_Genome2/data/6_longestTranscript/musNig/LongestTranscript.fasta  --out ~/Sorex_Genome2/analysis/1_busco/musNig/ -l mammalia_odb10 -m prot  -f

#note, longest transcript remove duplicated but also prunes transcripts that are a part of the busco score. so full busco run all
busco --in ~/Sorex_Genome2/data/5_TOGA/mSorAra2/query_prot.fasta  --out ~/Sorex_Genome2/analysis/1_busco/mSorAra2/ -l mammalia_odb10 -m prot  -f
#note, if first sequence does not start with M, busco gets confused thinking it is not a protein file, just swap the first sequence
