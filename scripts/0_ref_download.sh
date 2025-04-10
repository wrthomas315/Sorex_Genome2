#!/bin/bash

#Sorex files - used mSorAra2.pri.cur.20220707.fasta.gz from https://www.genomeark.org/vgp-curated-assembly/Sorex_araneus.html
#However, it is the same assembly deposited onto GenBank mSorAra2: GCA_027595985.1 and RefSeq GCF_027595985.1 and now easier to get instead of using AWS
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/595/985/GCF_027595985.1_mSorAra2.pri/GCF_027595985.1_mSorAra2.pri_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/595/985/GCF_027595985.1_mSorAra2.pri/GCF_027595985.1_mSorAra2.pri_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/595/985/GCF_027595985.1_mSorAra2.pri/GCF_027595985.1_mSorAra2.pri_rna.fna.gz
#
mv GCF_027595985.1_mSorAra2.pri_genomic.fna.gz mSorAra2.pri.cur.20220707.fasta.gz

#Human references
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
mv GCF_000001405.26_GRCh38_genomic.fna.gz hg38.fa.gz

#Mustela nigripes and Suncus etruscus for annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/355/385/GCF_022355385.1_MUSNIG.SB6536/GCF_022355385.1_MUSNIG.SB6536_genomic.fna.gz

#Eulipotyphla genomes for chromosome evolution and quast
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/115/265/GCA_004115265.3_mRhiFer1_v1.p/GCA_004115265.3_mRhiFer1_v1.p_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/126/475/GCA_004126475.3_mPhyDis1.pri.v3/GCA_004126475.3_mPhyDis1.pri.v3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/950/295/315/GCA_950295315.1_mEriEur2.1/GCA_950295315.1_mEriEur2.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/296/755/GCA_000296755.1_EriEur2.0/GCA_000296755.1_EriEur2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/181/395/GCA_000181395.1_ASM18139v1/GCA_000181395.1_ASM18139v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/139/225/GCA_024139225.1_mSunEtr1.pri.cur/GCA_024139225.1_mSunEtr1.pri.cur_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/363/575/GCA_004363575.1_SolPar_v1_BIUU/GCA_004363575.1_SolPar_v1_BIUU_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/898/055/GCA_014898055.3_MPIMG_talOcc4v2.1/GCA_014898055.3_MPIMG_talOcc4v2.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/455/555/GCA_019455555.1_Gpyr_1.0/GCA_019455555.1_Gpyr_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/260/355/GCA_000260355.1_ConCri1.0/GCA_000260355.1_ConCri1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/024/945/GCA_004024945.1_UroGra_v1_BIUU/GCA_004024945.1_UroGra_v1_BIUU_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/834/395/GCA_029834395.1_SorCin_2.0/GCA_029834395.1_SorCin_2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/021/461/705/GCA_021461705.1_Cryptotis_parva_assembly_1.0/GCA_021461705.1_Cryptotis_parva_assembly_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/024/925/GCA_004024925.1_ScaAqu_v1_BIUU/GCA_004024925.1_ScaAqu_v1_BIUU_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/027/635/GCA_004027635.1_CroInd_v1_BIUU/GCA_004027635.1_CroInd_v1_BIUU_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/565/675/GCA_028565675.1_ASM2856567v1/GCA_028565675.1_ASM2856567v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/901/085/GCA_002901085.1_ASM290108v1/GCA_002901085.1_ASM290108v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/181/275/GCA_000181275.2_SorAra2.0/GCA_000181275.2_SorAra2.0_genomic.fna.gz
mv *.gz ../data/0_refs/
