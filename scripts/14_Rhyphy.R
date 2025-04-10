# Goals
#STEP1: Look at BUSCOS for species used and plot
#STEP2: Import absrel results for foreground lineages and look at signifigance
#STEP3: Import absrel background selection
#STEP4: Subset them based on convergence and amount of background selection, can see what is shrew specific and Dehnel's related
#STEP5:. Plot Gene set enrichment for Figure 2 of shrew specific genes

#### Load in libraries needed
library(readr)
library(ggimage)
library(phyloseq)
library(ggtree)
library(TDbook)
library(tidyr)
library(tidyverse)
library(rphylopic)
library(ggmsa)
library(Biostrings)
library(UpSetR)
library(biomaRt)
library(fgsea)
library(ggstance)
library(stringr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(tibble)
library(AnnotationHub)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(gplots)

#STEP1: Look at BUSCOS for species used and plot
###Begin by looking at our buscos across the phylogeny
###Make graphs of BUSCOs and LT of species chosen
compGen_tree <- read.tree("~/Sorex_Genome2/data/9_trees/Alvarez_Carretero2022_4075sp.pruned.nh")
compGen_tree_images <- ggimage::phylopic_uid(compGen_tree$tip.label)
compGen_tree_images$order <- c(rep("Carnivora",8),rep("Artiodactyla",6),rep("Perissodactyla",1),rep("Pholidota",1),rep("Chiroptera",4),rep("Eulipotyphla",7),rep("Rodentia",5),rep("Lagomorpha",1),rep("Scadentia",1),rep("Primates",5))
compGen_tree_images$labeling <-  compGen_tree_images$name
compGen_tree_images$fore <-  c(rep("Back",1),rep("Fore",1),rep("Back",1),rep("Fore",1),rep("Back",17),rep("Fore",1),rep("Back",17))
compGen_tree_images$labeling
###one with species colored
species_colored <- ggtree(compGen_tree,branch.length="none") %<+% compGen_tree_images +
  geom_tiplab(aes(label=labeling,colour=order), offset=0, size =2)# + #xlim(NA, 1000)
###one with shrews and mustella labeled
fore_colored <- ggtree(compGen_tree) %<+% compGen_tree_images +
  geom_tiplab(aes(label=labeling,colour=fore), offset=0, size =2)+scale_color_manual(values=c("black","lightslateblue"))+theme_tree2()# + #xlim(NA, 1000) 
fore_colored
#####
CompGen_buscos <- read_table("~/Sorex_Genome2/analysis/1_busco/CompGen_buscos.txt")
CompGen_buscos <- as.data.frame(CompGen_buscos)
###Reshape the data to long format
cgp3 <- as.data.frame(cbind(CompGen_buscos$Species,CompGen_buscos$Miss,CompGen_buscos$Frag,CompGen_buscos$BUSCO_complete))
cgp3$V2 <- as.numeric(cgp3$V2)
cgp3$V3 <- as.numeric(cgp3$V3)
cgp3$V4 <- as.numeric(cgp3$V4)
busco_long <- reshape2::melt(cgp3, id.vars = "V1",
                             measure.vars = c("V2", "V3", "V4"),
                             variable.name = "category", value.name = "percent")

fore_colored +  geom_facet(panel = "Busco", data = busco_long, geom = geom_barh, mapping = aes(x = percent, fill = as.factor(category)),stat='identity' )+scale_fill_manual(values = c("#BBDE28FF","#EBE51AFF","#2F6B8EFF"))+scale_y_continuous()+theme_tree2()
#####





#STEP2: Subset absrel results for foreground branches
###Positive selection###
###Get list from hyphy explore and hyphy foreground
abs_for <- read_table("~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/output.txt", 
                      col_names = FALSE)
#import transcript to gene and add to column to make a key
LT_output_4hp <- read_table("~/Sorex_Genome2/analysis/3_hyphy/1_abs_fore/LT_output_4hp.txt", 
                            col_names = FALSE)
abs_for$X5 <- NA
for (i in seq_len(nrow(abs_for))) {
  # Find the matching row
  match_row <- match(abs_for[i, "X1"], LT_output_4hp$X3)
  
  # If a match found, assign the corresponding value to X5
  if (!is.na(match_row)) {
    abs_for[i, "X5"] <- LT_output_4hp[match_row, "X1"]
  }
}
#use gene ID to get ensemble gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror = "useast")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
abs_for_vec <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = abs_for$X5, 
                     mart = ensembl)
abs_for$X6 <- NA
# Iterate over rows in absrel results
for (i in seq_len(nrow(abs_for))) {
  # Find the matching row in ensembl data frame
  match_row <- match(abs_for[i, "X5"], abs_for_vec$ensembl_gene_id)
  
  # If a match found, assign the corresponding value to X6
  if (!is.na(match_row)) {
    abs_for[i, "X6"] <- abs_for_vec[match_row, "external_gene_name"]
  }
}
#Remove NAs or genes without positive selection for each species and make vector
ara.sub <-subset(abs_for, !is.na(X3) & X2 == "sorex_araneus")
etr.sub <-subset(abs_for, !is.na(X3) & X2 == "suncus_etruscus")
fur.sub <-subset(abs_for, !is.na(X3) & X2 == "mustela_putorius_furo")
erm.sub <-subset(abs_for, !is.na(X3) & X2 == "mustela_erminea")
#Can see how many genes overlap between our four foreground species
gene_listUPSET <- list(Mustela_erminea = erm.sub$X1,Sorex_araneus = ara.sub$X1,Mustela_putorius_furo = fur.sub$X1,Suncus_etruscus = etr.sub$X1)
upset_gene_listUPSET <- fromList(gene_listUPSET)
upset(upset_gene_listUPSET,nsets = 4,order.by = "freq",)
#How many are significant in each species without correction (testing four foreground branches across all genes)
ara.sub2 <-subset(abs_for, !is.na(X3) & X2 == "sorex_araneus" & X3 <= .05)
etr.sub2 <-subset(abs_for, !is.na(X3) & X2 == "suncus_etruscus" & X3 <= .05)
fur.sub2 <-subset(abs_for, !is.na(X3) & X2 == "mustela_putorius_furo" & X3 <= .05)
erm.sub2 <-subset(abs_for, !is.na(X3) & X2 == "mustela_erminea" & X3 <= .05)
gene_listUPSET2 <- list(Mustela_erminea = erm.sub2$X1,Sorex_araneus = ara.sub2$X1,Mustela_putorius_furo = fur.sub2$X1,Suncus_etruscus = etr.sub2$X1)
upset_gene_listUPSET2 <- fromList(gene_listUPSET2)
upset(upset_gene_listUPSET2,nsets = 4,order.by = "freq",)
Reduce(intersect, gene_listUPSET2)
#How many are significant in each species with only one correction (testing four foreground branches, not correcting across all genes tested)
ara.sub3 <-subset(abs_for, !is.na(X3) & X2 == "sorex_araneus" & X4 <= .05)
etr.sub3 <-subset(abs_for, !is.na(X3) & X2 == "suncus_etruscus" & X4 <= .05)
fur.sub3 <-subset(abs_for, !is.na(X3) & X2 == "mustela_putorius_furo" & X4 <= .05)
erm.sub3 <-subset(abs_for, !is.na(X3) & X2 == "mustela_erminea" & X4 <= .05)
gene_listUPSET3 <- list(Mustela_erminea = erm.sub3$X1,Sorex_araneus = ara.sub3$X1,Mustela_putorius_furo = fur.sub3$X1,Suncus_etruscus = etr.sub3$X1)
upset_gene_listUPSET3 <- fromList(gene_listUPSET3)
upset(upset_gene_listUPSET3,nsets = 4,order.by = "freq",)
#Finally, wow many are significant in each species with both corrections (testing four foreground branches and correcting across all genes tested)
ara.sub4 <-subset(abs_for, !is.na(X3) & X2 == "sorex_araneus")
etr.sub4 <-subset(abs_for, !is.na(X3) & X2 == "suncus_etruscus")
fur.sub4 <-subset(abs_for, !is.na(X3) & X2 == "mustela_putorius_furo")
erm.sub4 <-subset(abs_for, !is.na(X3) & X2 == "mustela_erminea")
ara.sub4$p_adjusted <- p.adjust(ara.sub4$X4, method = "BH")
etr.sub4$p_adjusted <- p.adjust(etr.sub4$X4, method = "BH")
fur.sub4$p_adjusted <- p.adjust(fur.sub4$X4, method = "BH")
erm.sub4$p_adjusted <- p.adjust(erm.sub4$X4, method = "BH")
ara.sub5 <-subset(ara.sub4, p_adjusted <= .05)
etr.sub5 <-subset(etr.sub4,p_adjusted <= .05)
fur.sub5 <-subset(fur.sub4, p_adjusted <= .05)
erm.sub5 <-subset(erm.sub4, p_adjusted <= .05)
gene_listUPSET5 <- list(Mustela_erminea = erm.sub5$X1,Sorex_araneus = ara.sub5$X1,Mustela_putorius_furo = fur.sub5$X1,Suncus_etruscus = etr.sub5$X1)
upset_gene_listUPSET5 <- fromList(gene_listUPSET5)
upset(upset_gene_listUPSET5,nsets = 4,order.by = "freq",)
Reduce(intersect, gene_listUPSET5)





###STEP 3: Look at background selection to see what is shrew or Dehnel's specific
###Now adding in background selection to make species/Dehnel's specific
#excluding other branches no nodes, need to add correction for all genes tested per species
abs_exp <- read_table("~/Sorex_Genome2/analysis/3_hyphy/2_abs_dropout/output_noforenode.txt", 
                      col_names = FALSE)
#homo_sapiens
homSap.sub <-subset(abs_exp, !is.na(X3) & X2 == "homo_sapiens")
homSap.sub$p_adjusted <- p.adjust(homSap.sub$X4, method = "BH")
homSap.sub2 <-subset(homSap.sub, p_adjusted <= .05)
#phyllostomus_discolor
phyDis.sub <-subset(abs_exp, !is.na(X3) & X2 == "phyllostomus_discolor")
phyDis.sub$p_adjusted <- p.adjust(phyDis.sub$X4, method = "BH")
phyDis.sub2 <-subset(phyDis.sub, p_adjusted <= .05)
#acomys_russatus
acoRus.sub <-subset(abs_exp, !is.na(X3) & X2 == "acomys_russatus")
acoRus.sub$p_adjusted <- p.adjust(acoRus.sub$X4, method = "BH")
acoRus.sub2 <-subset(acoRus.sub, p_adjusted <= .05)
#talpa_occidentalis
talOcc.sub <-subset(abs_exp, !is.na(X3) & X2 == "talpa_occidentalis")
talOcc.sub$p_adjusted <- p.adjust(talOcc.sub$X4, method = "BH")
talOcc.sub2 <-subset(talOcc.sub, p_adjusted <= .05)
#tupaia_chinensis
tupChi.sub <-subset(abs_exp, !is.na(X3) & X2 == "tupaia_chinensis")
tupChi.sub$p_adjusted <- p.adjust(tupChi.sub$X4, method = "BH")
tupChi.sub2 <-subset(tupChi.sub, p_adjusted <= .05)
#microcebus_murinus
micMur.sub <-subset(abs_exp, !is.na(X3) & X2 == "microcebus_murinus")
micMur.sub$p_adjusted <- p.adjust(micMur.sub$X4, method = "BH")
micMur.sub2 <-subset(micMur.sub, p_adjusted <= .05)
#papio_anubis
papAnu.sub <-subset(abs_exp, !is.na(X3) & X2 == "papio_anubis")
papAnu.sub$p_adjusted <- p.adjust(papAnu.sub$X4, method = "BH")
papAnu.sub2 <-subset(papAnu.sub, p_adjusted <= .05)
#hylobates_moloch
hylMol.sub <-subset(abs_exp, !is.na(X3) & X2 == "hylobates_moloch")
hylMol.sub$p_adjusted <- p.adjust(hylMol.sub$X4, method = "BH")
hylMol.sub2 <-subset(hylMol.sub, p_adjusted <= .05)
#pan_troglodytes
panTro.sub <-subset(abs_exp, !is.na(X3) & X2 == "pan_troglodytes")
panTro.sub$p_adjusted <- p.adjust(panTro.sub$X4, method = "BH")
panTro.sub2 <-subset(panTro.sub, p_adjusted <= .05)
#macaca_mulatta
macMul.sub <-subset(abs_exp, !is.na(X3) & X2 == "macaca_mulatta")
macMul.sub$p_adjusted <- p.adjust(macMul.sub$X4, method = "BH")
macMul.sub2 <-subset(macMul.sub, p_adjusted <= .05)
#peromyscus_leucopus
perLeu.sub <-subset(abs_exp, !is.na(X3) & X2 == "peromyscus_leucopus")
perLeu.sub$p_adjusted <- p.adjust(perLeu.sub$X4, method = "BH")
perLeu.sub2 <-subset(perLeu.sub, p_adjusted <= .05)
#rattus_norvegicus
ratNor.sub <-subset(abs_exp, !is.na(X3) & X2 == "rattus_norvegicus")
ratNor.sub$p_adjusted <- p.adjust(ratNor.sub$X4, method = "BH")
ratNor.sub2 <-subset(ratNor.sub, p_adjusted <= .05)
#sciurus_vulgaris
sciVul.sub <-subset(abs_exp, !is.na(X3) & X2 == "sciurus_vulgaris")
sciVul.sub$p_adjusted <- p.adjust(sciVul.sub$X4, method = "BH")
sciVul.sub2 <-subset(sciVul.sub, p_adjusted <= .05)
#oryctolagus_cuniculus
oryCun.sub <-subset(abs_exp, !is.na(X3) & X2 == "oryctolagus_cuniculus")
oryCun.sub$p_adjusted <- p.adjust(oryCun.sub$X4, method = "BH")
oryCun.sub2 <-subset(oryCun.sub, p_adjusted <= .05)
#pipistrellus_pipistrellus
pipPip.sub <-subset(abs_exp, !is.na(X3) & X2 == "pipistrellus_pipistrellus")
pipPip.sub$p_adjusted <- p.adjust(pipPip.sub$X4, method = "BH")
pipPip.sub2 <-subset(pipPip.sub, p_adjusted <= .05)
#myotis_myotis
myoMyo.sub <-subset(abs_exp, !is.na(X3) & X2 == "myotis_myotis")
myoMyo.sub$p_adjusted <- p.adjust(myoMyo.sub$X4, method = "BH")
myoMyo.sub2 <-subset(myoMyo.sub, p_adjusted <= .05)
#rhinolophus_ferrumequinum
rhiFer.sub <-subset(abs_exp, !is.na(X3) & X2 == "rhinolophus_ferrumequinum")
rhiFer.sub$p_adjusted <- p.adjust(rhiFer.sub$X4, method = "BH")
rhiFer.sub2 <-subset(rhiFer.sub, p_adjusted <= .05)
#camelus_ferus
camFer.sub <-subset(abs_exp, !is.na(X3) & X2 == "camelus_ferus")
camFer.sub$p_adjusted <- p.adjust(camFer.sub$X4, method = "BH")
camFer.sub2 <-subset(camFer.sub, p_adjusted <= .05)
#bos_gaurus
bosGau.sub <-subset(abs_exp, !is.na(X3) & X2 == "bos_gaurus")
bosGau.sub$p_adjusted <- p.adjust(bosGau.sub$X4, method = "BH")
bosGau.sub2 <-subset(bosGau.sub, p_adjusted <= .05)
#ovis_orientalis
oviOri.sub <-subset(abs_exp, !is.na(X3) & X2 == "ovis_orientalis")
oviOri.sub$p_adjusted <- p.adjust(oviOri.sub$X4, method = "BH")
oviOri.sub2 <-subset(oviOri.sub, p_adjusted <= .05)
#phocoena_sinus
phoSin.sub <-subset(abs_exp, !is.na(X3) & X2 == "phocoena_sinus")
phoSin.sub$p_adjusted <- p.adjust(phoSin.sub$X4, method = "BH")
phoSin.sub2 <-subset(phoSin.sub, p_adjusted <= .05)
#sus_scrofa
susScr.sub <-subset(abs_exp, !is.na(X3) & X2 == "sus_scrofa")
susScr.sub$p_adjusted <- p.adjust(susScr.sub$X4, method = "BH")
susScr.sub2 <-subset(susScr.sub, p_adjusted <= .05)
#ceratotherium_simum_cottoni
cerSim.sub <-subset(abs_exp, !is.na(X3) & X2 == "ceratotherium_simum_cottoni")
cerSim.sub$p_adjusted <- p.adjust(cerSim.sub$X4, method = "BH")
cerSim.sub2 <-subset(cerSim.sub, p_adjusted <= .05)
#manis_javanica
manJav.sub <-subset(abs_exp, !is.na(X3) & X2 == "manis_javanica")
manJav.sub$p_adjusted <- p.adjust(manJav.sub$X4, method = "BH")
manJav.sub2 <-subset(manJav.sub, p_adjusted <= .05)
#vulpes_lagopus
vulLag.sub <-subset(abs_exp, !is.na(X3) & X2 == "vulpes_lagopus")
vulLag.sub$p_adjusted <- p.adjust(vulLag.sub$X4, method = "BH")
vulLag.sub2 <-subset(vulLag.sub, p_adjusted <= .05)
#lynx_canadensis
lynCan.sub <-subset(abs_exp, !is.na(X3) & X2 == "lynx_canadensis")
lynCan.sub$p_adjusted <- p.adjust(lynCan.sub$X4, method = "BH")
lynCan.sub2 <-subset(lynCan.sub, p_adjusted <= .05)
#ailuropoda_melanoleuca
ailMel.sub <-subset(abs_exp, !is.na(X3) & X2 == "ailuropoda_melanoleuca")
ailMel.sub$p_adjusted <- p.adjust(ailMel.sub$X4, method = "BH")
ailMel.sub2 <-subset(ailMel.sub, p_adjusted <= .05)
#mustela_nigripes
musNig.sub <-subset(abs_exp, !is.na(X3) & X2 == "mustela_nigripes")
musNig.sub$p_adjusted <- p.adjust(musNig.sub$X4, method = "BH")
musNig.sub2 <-subset(musNig.sub, p_adjusted <= .05)
#lutra_lutra
lutLut.sub <-subset(abs_exp, !is.na(X3) & X2 == "lutra_lutra")
lutLut.sub$p_adjusted <- p.adjust(lutLut.sub$X4, method = "BH")
lutLut.sub2 <-subset(lutLut.sub, p_adjusted <= .05)
#martes_zibellina
marZib.sub <-subset(abs_exp, !is.na(X3) & X2 == "martes_zibellina")
marZib.sub$p_adjusted <- p.adjust(marZib.sub$X4, method = "BH")
marZib.sub2 <-subset(marZib.sub, p_adjusted <= .05)
#condylura_cristata
conCri.sub <-subset(abs_exp, !is.na(X3) & X2 == "condylura_cristata")
conCri.sub$p_adjusted <- p.adjust(conCri.sub$X4, method = "BH")
conCri.sub2 <-subset(conCri.sub, p_adjusted <= .05)
#erinaceus_europaeus
eriEur.sub <-subset(abs_exp, !is.na(X3) & X2 == "erinaceus_europaeus")
eriEur.sub$p_adjusted <- p.adjust(eriEur.sub$X4, method = "BH")
eriEur.sub2 <-subset(eriEur.sub, p_adjusted <= .05)
#scalopus_aquaticus
scaAqu.sub <-subset(abs_exp, !is.na(X3) & X2 == "scalopus_aquaticus")
scaAqu.sub$p_adjusted <- p.adjust(scaAqu.sub$X4, method = "BH")
scaAqu.sub2 <-subset(scaAqu.sub, p_adjusted <= .05)
#solenodon_paradoxus
solPar.sub <-subset(abs_exp, !is.na(X3) & X2 == "solenodon_paradoxus")
solPar.sub$p_adjusted <- p.adjust(solPar.sub$X4, method = "BH")
solPar.sub2 <-subset(solPar.sub, p_adjusted <= .05)
#mus_musculus
musMus.sub <-subset(abs_exp, !is.na(X3) & X2 == "mus_musculus")
musMus.sub$p_adjusted <- p.adjust(musMus.sub$X4, method = "BH")
musMus.sub2 <-subset(musMus.sub, p_adjusted <= .05)
#cervus_elaphus
cerEla.sub <-subset(abs_exp, !is.na(X3) & X2 == "cervus_elaphus")
cerEla.sub$p_adjusted <- p.adjust(cerEla.sub$X4, method = "BH")
cerEla.sub2 <-subset(cerEla.sub, p_adjusted <= .05)
#
combined_vector <- rbind(cerEla.sub2,musMus.sub2,solPar.sub2,scaAqu.sub2,eriEur.sub2,conCri.sub2,marZib.sub2,homSap.sub2,phyDis.sub2,acoRus.sub2,talOcc.sub2,tupChi.sub2,micMur.sub2,papAnu.sub2,hylMol.sub2,panTro.sub2,macMul.sub2,perLeu.sub2,ratNor.sub2,sciVul.sub2,oryCun.sub2,pipPip.sub2,myoMyo.sub2,rhiFer.sub2,camFer.sub2,bosGau.sub2,oviOri.sub2,phoSin.sub2,susScr.sub2,cerSim.sub2,vulLag.sub2,manJav.sub2,lynCan.sub2,musNig.sub2,lutLut.sub2,ailMel.sub2)
unique_vector <- unique(combined_vector$X1)
gene_listUPSET6 <- list(Mustela_erminea = erm.sub5$X1,Sorex_araneus = ara.sub5$X1,Mustela_putorius_furo = fur.sub5$X1,Suncus_etruscus = etr.sub5$X1,Other =unique_vector)
upset_gene_listUPSET6 <- fromList(gene_listUPSET6)
upset(upset_gene_listUPSET6,nsets = 5,order.by = "freq",)
Reduce(intersect, gene_listUPSET6)
subset(combined_vector, !is.na(X3) & X1 == "ENST00000348564")
#gene found in mustelids and shrew but not other column
setdiff(Reduce(intersect, list(gene_listUPSET6$Mustela_erminea, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_putorius_furo)),gene_listUPSET6$Other)
#Less strict convergence with no background selection
enst_counts <- table(combined_vector$X1)
result_gt5 <- names(enst_counts[enst_counts > 5])
gene_listUPSET7 <- list(Mustela_erminea = erm.sub5$X1,Sorex_araneus = ara.sub5$X1,Mustela_putorius_furo = fur.sub5$X1,Suncus_etruscus = etr.sub5$X1,OtherGT5 =result_gt5)
upset_gene_listUPSET7 <- fromList(gene_listUPSET7)
upset(upset_gene_listUPSET7,nsets = 5,order.by = "freq",)
setdiff(Reduce(intersect, list(gene_listUPSET6$Suncus_etruscus, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_erminea)),gene_listUPSET7$OtherGT5)
###okay we see the numbers now. lets take a look at what things/bins contain
#what is the distribution
ggplot(ara.sub4, aes(x=p_adjusted)) + geom_density()+xlim(-1,15)
plot(ecdf(ara.sub4$p_adjusted))





###STEP 4: Subset them based on convergence and amount of background selection and write ouputs
#676 shrew (ara.sub5)
write.table(ara.sub5, file='~/Sorex_Genome2/analysis/3_hyphy/shrewtotal.txt', quote=FALSE, sep='\t')
#231 shrew specifc
shrewspecific <-subset(ara.sub5, !(X1 %in% combined_vector$X1) & !(X1 %in% erm.sub5$X1)& !(X1 %in% etr.sub5$X1)& !(X1 %in% fur.sub5$X1))
write.table(shrewspecific, file='~/Sorex_Genome2/analysis/3_hyphy/shrewspecific.txt', quote=FALSE, sep='\t')
#23 convergent regardless of background, then with gene names
CONV23 <-unique(c(Reduce(intersect, list(gene_listUPSET6$Mustela_erminea, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_putorius_furo)),Reduce(intersect, list(gene_listUPSET6$Suncus_etruscus, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_putorius_furo)),Reduce(intersect, list(gene_listUPSET6$Mustela_erminea, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Suncus_etruscus))))
GNCONV23<-unique(c(Reduce(intersect, list(erm.sub5$X6, ara.sub5$X6,fur.sub5$X6)),Reduce(intersect, list(etr.sub5$X6, ara.sub5$X6,fur.sub5$X6)),Reduce(intersect, list(erm.sub5$X6, ara.sub5$X6,etr.sub5$X6))))
#10 convergent loose background, then with gene names
looseCONV8 <-unique(c(setdiff(Reduce(intersect, list(gene_listUPSET6$Suncus_etruscus, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_erminea)),gene_listUPSET7$OtherGT5),setdiff(Reduce(intersect, list(gene_listUPSET6$Mustela_putorius_furo, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_erminea)),gene_listUPSET7$OtherGT5),setdiff(Reduce(intersect, list(gene_listUPSET6$Suncus_etruscus, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_putorius_furo)),gene_listUPSET7$OtherGT5)))
combined_vector$X5 <- NA
for (i in seq_len(nrow(combined_vector))) {
  # Find the matching row
  match_row <- match(combined_vector[i, "X1"], LT_output_4hp$X3)
  
  # If a match found, assign the corresponding value to X5
  if (!is.na(match_row)) {
    combined_vector[i, "X5"] <- LT_output_4hp[match_row, "X1"]
  }
}
#use gene ID to get ensemble gene name
combined_vector_vec <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                             filters = 'ensembl_gene_id',
                             values = combined_vector$X5, 
                             mart = ensembl)
combined_vector$X6 <- NA
# Iterate over rows in absrel results
for (i in seq_len(nrow(combined_vector))) {
  # Find the matching row in ensembl data frame
  match_row <- match(combined_vector[i, "X5"], combined_vector_vec$ensembl_gene_id)
  
  # If a match found, assign the corresponding value to X6
  if (!is.na(match_row)) {
    combined_vector[i, "X6"] <- combined_vector_vec[match_row, "external_gene_name"]
  }
}
GMenst_counts <- table(combined_vector$X6)
GMresult_gt5 <- names(GMenst_counts[GMenst_counts > 5])
GNlooseCONV8 <-unique(c(setdiff(Reduce(intersect, list(etr.sub5$X6, ara.sub5$X6,erm.sub5$X6)),GMresult_gt5),setdiff(Reduce(intersect, list(fur.sub5$X6, ara.sub5$X6,erm.sub5$X6)),GMresult_gt5),setdiff(Reduce(intersect, list(etr.sub5$X6, ara.sub5$X6,fur.sub5$X6)),GMresult_gt5)))
setdiff(setdiff(Reduce(intersect, list(erm.sub5$X6, ara.sub5$X6,fur.sub5$X6)),GMresult_gt5),etr.sub5$X6)
#1 (PCDHA6) convergent strict background
setdiff(Reduce(intersect, list(gene_listUPSET6$Mustela_erminea, gene_listUPSET6$Sorex_araneus,gene_listUPSET6$Mustela_putorius_furo)),gene_listUPSET6$Other)
#1 (SLC2A11) convergent in all 4 but found in background 6 times
Reduce(intersect, gene_listUPSET5)
#Final upset for supplement
gene_listUPSET8 <- list(Mustela_erminea = erm.sub5$X1,Sorex_araneus = ara.sub5$X1,Mustela_putorius_furo = fur.sub5$X1,Suncus_etruscus = etr.sub5$X1,Loose =result_gt5, Strict=unique_vector)
upset_gene_listUPSET8 <- fromList(gene_listUPSET8)
upset(upset_gene_listUPSET8,nsets = 6,order.by = "freq",)
#total list for meme
total_vector <- rbind(combined_vector,ara.sub5,etr.sub5,fur.sub5,erm.sub5)
write.table(total_vector, file='~/Sorex_Genome2/analysis/3_hyphy/hyphytotal.txt', quote=FALSE, sep='\t')
#
length(GMresult_gt5)
length(GMenst_counts)
GMenst_counts
#
length(result_gt5)
result_gt5
length(unique_vector)
unique_vector


###STETP 5: Pathwya enrichment of shrew specific psgs
library(viridis)

#########
#KEGG for david clusters, total shrew
total_workaable <- read_delim("~/Sorex_Genome2/analysis/3_hyphy/3_KEGG/KEGG_shrewtotal_workable.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
total_workaable$Term <- factor(total_workaable$Term, levels = total_workaable$Term[order(-total_workaable$PValue)])
fanconi_png <- ggplot(total_workaable , aes(Term,-log(PValue),color = Count))+
  geom_point(size=5)+
  #scale_color_viridis(option = "C")+
  coord_flip()+
  theme_bw()+theme(axis.text.y = element_blank(),  # Remove y-axis text labels
                   axis.title.y = element_blank())
ggsave("~/Sorex_Genome2/analysis/3_hyphy/3_KEGG/GSEA.png", fanconi_png,width = 4, height = 5, dpi =300,)
ggplot(total_workaable , aes(Term,-log(PValue),color = Count))+
  geom_point(size=5)+
  #scale_color_viridis(option = "C")+
  coord_flip()+
  theme_bw()


