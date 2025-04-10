## *Sorex araneus* Genome Paper
Scripts and code to reproduce RNAseq assembly, annotation, and analyses of the common shrew, *Sorex araneus*, reference genome.

### Goals and strategy

The objective of this project is to infer the evolutionary processes and functional mechanisms underlying the shrews unique trait set, which includes seasonal size plasticity, extensive intraspecies chromosomal rearrangements, an extremely short lifespan, and 
extraordinarily high metabolic rates. To conduct this research, we first assembled a chromosomal reference genome for *Sorex araneus* with the Vertebrate Genomes Project and compared chromosomes to species within the order Eulipotyphla. Then, using novel 
annotations and alignments from TOGA, we tested genes for positive selection to identify candidate genes with aBSREL and MEME. Finally, we compared these candidate genes to those identified in previous and current functional analyses.

Note1: There are various spots you could start at with this analysis using this github. You can download the *S. araneus* reference genome, mask it, align it the human genome, annotate it, align it to the other mammal species genes, clean the 
alignments, and then run evolutionary analyses. Alternatively, you could grab alignments supplied here and continue to hyphy analyses. Same thing goes for the RNA seq, you could start with cleaning and alignment, or use presupplied files for DESeq2 
differential expression analyses. Your choice!

Note2: Some of the most difficult parts of reproducing these experiments is having the proper environments, cluster space, and cluster memory allocations to run things, such as TOGA. I will advise on what programs were in my local 
environment throughout the process. However, code to download each program will be found on their respective githubs.

### Data, Code, and Analysis

#### Assembly

The new *S.araneus* reference genome was assembled with the Vertebrate Genomes Project workflow (v2.0), which can be found here https://github.com/VGP/vgp-assembly.

#### Download assembly, transcript files, and other references. Set up environments.

Get *S. araneus* assembly/ncbi transcriptome, human reference to compare it to, as well as other genomes for chromosome number comparisons and *Mustela nigripes* annotation.
Note: Also placed in folder already is human reference file used for TOGA annotations.

```
bash 0_ref_download.sh
bash 1_buildenvs.sh
```

#### Masking genomes
Note: Make sure to go into scripts in scripts_bin to change cluster settings, user name in 2_mask.sh. Masked genomes in ~/Sorex_Genome2/data/1_mask
```
bash 2_mask.sh
```

#### Genome alignments
Note: Would run each species separately and need to configure to individual cluster
```
bash 3_alignment.sh
```


#### Cleaning genome alignments
Clean genome alignments in ~/Sorex_Genome2/data/4_clean_genome_align
```
bash 4_cleaning_align.sh
```

#### Annotation with TOGA
Note: Make sure you have Toga installed, also activate Toga environment.
```
bash 5_TOGA.sh
```

#### Get your longest transcripts
```
bash 6_longest_transcript.sh          
```

#### Test quality of TOGA with BUSCO
Primary buscos in ~/Sorex_Genome2/data/7_busco with complete data table here
```
bash 7_busco.sh          
```

#### Gene alignment with MACSE exon by exon
Note: only including taper cleaned alignments below to reduce meta data memory storage
```
bash 8_exonbyexon.sh            
```

#### Clean alignments with TAPER
Gene alignments in zip file ~/Sorex_Genome2/data/9_taper_align
```
bash 9_taper.sh            
```

#### Prepare primary files for use in hyphy
Prune unrooted trees for each individual alignment.
```
python 10_treeprune.py            
```
Modify taper alignments to set foreground lineage
```
11_hyphy_prep.sh            
```


#### Prepare aBSREL jobList for foreground and background dropout runs then run in parallel on own
```
bash 12_absrel.sh          
```
Convert the output from jsons
```
bash 13
```
And then run in R to test for significance
```
R 14_Rhyphy.R   
```
You will receive:
All Sorex araneus PSGs ~/Sorex_Genome2/analysis/3_hyphy/shrewtotal.txt
Sorex araneus specific ~/Sorex_Genome2/analysis/3_hyphy/PSGs/shrewspecific.txt
Genes and species to find PSH sites in meme ~/Sorex_Genome2/analysis/3_hyphy/hyphytotal.txt
Overlapping genes within script
Note: KEGG ran on web browser, results in ~/Sorex_Genome2/analysis/3_hyphy/3_KEGG

#### MEME
Create trees for meme, labeling foreground branches with species identified to have gene under positive selection from absrel. Then run meme for those species and genes from jobList.
```
bash 15_meme.sh            
```
Then sort the jsons to identify which sites are significant and convergent between species tested.
```
bash 16_meme_output.sh            
```
Will end up with MEMEoutput.txt for all sites for species tested and convergencesiteShrew.txt for convergent sites.

#### Transcriptomics
Need to clean sequencing with fastp, pseudo-align to reference transcriptome with Kalisto
```
bash 17_fastp_kallisto.sh           
```
Differential expression of cortex and hippocampus between seasons and comparisons to positively selected genes.
```
R 18_Cortex_Hippocampus.R          
```
Then a couple of additional scripts from the results that will make figures for interesting genes. Best done in R studio.
```
R 19_PSG_upset_fig.R
R 20_VEGF_figure.R
R 21_PCDHA6_fig.R
```
Note 19 also runs permutation test for overlap of convergent genes.


#### Quast
```
bash 22_quast.sh
```


#### Chromosome evolution
Have genome alignments from Steps 2-4.
```
23_deschrambler.txt
```
and visualize (best done in Rstudio)
```
R 24_deschram_vis.R
```

