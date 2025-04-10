#!/bin/bash

#get TAPER
wget 
https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.4-linux-x86_64.tar.gz ./bin
tar zxvf julia-1.9.4-linux-x86_64.tar.gz

#parallelize and run taper on each alignment
module load gnu-parallel/6.0
parallel --verbose -j 40 --link "julia ./bin/julia-1.9.4/bin/correction_multi.jl -m N -a N {1} > ~/ShrewProjects/Sorex_Genome2/data/8_taper/{1}.taper.fa" ::: ~/ShrewProjects/Sorex_Genome2/data/7_alignments/*.macse.aligned.fa

