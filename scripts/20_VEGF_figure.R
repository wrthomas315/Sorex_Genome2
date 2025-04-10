###Making figure for comparing the hypothalamus data to sherew data
library(stringr)
library(readr)
library(assertr)
library(pheatmap)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(tibble)
library(tidyr)
library(AnnotationHub)
library(pheatmap)
library(RColorBrewer)
library(ggtree)
library(ggimage)
library(phyloseq)
library(ggtree)
library(TDbook)
library(tidyr)
library(rphylopic)
library(dplyr)
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



###LOAD in transcripts and make heatmap
####loading in data make a fun little graph
adult_nopap_hypoth <- read.table("~/Sorex_Genome2/data/000_miscFiles/Hypothal_data_2024/adult_nopap_hypothALV", quote="\"", comment.char="")
adult_nopap_hypot <- adult_nopap_hypoth[,-1]
rownames(adult_nopap_hypot) <- adult_nopap_hypoth[,1]
maybe <- t(adult_nopap_hypot)
maybe <- as.data.frame(maybe)
SPECIES <- c("capHir","oviAri","bosTur","susScr","sorAra","cavApe","hetGla","fukMec","perMan","micOch","musMus","ratNor","speTri","macMul","homSap","panTro")
SPECIES2 <- factor(c(rep("capHir", 3),rep("oviAri",7),rep("bosTur",6),rep("susScr",8),rep("sorAra",5),rep("cavApe",8),rep("hetGla",6),rep("fukMec",8),rep("perMan",8),rep("micOch",8),rep("musMus",4),rep("ratNor",7),rep("speTri",5),rep("macMul",7),rep("homSap",4),rep("panTro",3)))
SPECIES2 <- factor(SPECIES2, levels = SPECIES)
eve4graph <- cbind(maybe,as.data.frame(SPECIES2))
#and can take a loot at any gene of interest
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000009872)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()

#intersecting genes
#VEGFA ENSBTAG00000005339
#SPHK2 ENSBTAG00000009872
#MLST8 ENSBTAG00000009928
#GRIP2 ENSBTAG00000011896
#FANCI ENSBTAG00000009097
#KIAA0930 ENSBTAG00000006776
#LOXL1 ENSBTAG00000009086
#PIGR ENSBTAG00000019798
#GALNT12 ENSBTAG00000009037
#POLDIP3ENSBTAG00000016511
#NXPH3 ENSBTAG00000015293

###Heatmap of important gebes
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000005339,eve4graph$ENSBTAG00000009872,eve4graph$ENSBTAG00000009928,eve4graph$ENSBTAG00000011896,eve4graph$ENSBTAG00000009097,eve4graph$ENSBTAG00000006776,eve4graph$ENSBTAG00000009086,eve4graph$ENSBTAG00000019798,eve4graph$ENSBTAG00000009037,eve4graph$ENSBTAG00000016511,eve4graph$ENSBTAG00000015293))
tree <- read.tree("~/Sorex_Genome2/data/000_miscFiles/Hypothal_data_2024/hypot_ALVAREZ_TreePruned4ggtree.nh")
tree2 <- read.tree("~/Sorex_Genome2/data/000_miscFiles/Hypothal_data_2024/hypot_ALVAREZ_TreePruned.nh")

tree_key <- data.frame(key1=tree$tip.label,
                       key2=tree2$tip.label)
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
expDat$V2 <- as.numeric(expDat$V2)
expDat$V3 <- as.numeric(expDat$V3)
expDat$V4 <- as.numeric(expDat$V4)
expDat$V5 <- as.numeric(expDat$V5)
expDat$V6 <- as.numeric(expDat$V6)
expDat$V7 <- as.numeric(expDat$V7)
expDat$V8 <- as.numeric(expDat$V8)
expDat$V9 <- as.numeric(expDat$V9)
expDat$V10 <- as.numeric(expDat$V10)
expDat$V11 <- as.numeric(expDat$V11)
expDat$V12 <- as.numeric(expDat$V11)
mean_expDat1 <-aggregate(.~V1, data=expDat, mean)
hyp_heat <- as.data.frame(mean_expDat[,2:12])
hyp_heat2 <- sapply(hyp_heat, function(hyp_heat) (hyp_heat-mean(hyp_heat))/sd(hyp_heat))
rownames(hyp_heat2)<-mean_expDat1[,1]
colnames(hyp_heat2)<-c("VEGFA","SPHK2","MLST8","GRIP2","FANCI","KIAA0930","LOXL1","PIGR","GALNT12","POLDIP3","NXPH3")
hyp_heat2<-hyp_heat2[-c(2, 3, 4,5,8), ]
rownames(hyp_heat2) <- c("bos_gaurus","homo_sapiens","macaca_mulatta","mus_musculus","ovis_orientalis","pan_troglodytes","peromyscus_leucopus","rattus_norvegicus","sorex_araneus","sciurus_vulgaris","sus_scrofa")
hyp_heat2





#LOAD in desired tree
###going with VEGFA as has all 40 species 
VEGFjan2025Tree <- read.tree("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/VEGF_figure/tree/11spec_tree.nh")
#phylopic depreciated but still run to create needed frame for fore
VEGFjan2025Tree_images <- ggimage::phylopic_uid(VEGFjan2025Tree$tip.label)
VEGFjan2025Tree_images$labeling <-  VEGFjan2025Tree_images$name
VEGFjan2025Tree_images$labeling
VEGFjan2025Tree_images$fore <- c(rep("Background",3),rep("Foreground",1),rep("Background",7))
#check to make sure foreground and background are w/ correct species
VEGFjan2025Tree_images
VEGFjan2025Tree_phy <-ggtree(VEGFjan2025Tree,branch.length="none") %<+% VEGFjan2025Tree_images + 
  geom_tippoint(aes(label=labeling,color=fore),size=3, alpha=.5)+
  geom_tiplab(aes(label=labeling))+
  scale_color_manual(values=c("black","lightslateblue"))
VEGFjan2025Tree_phy


###LOAD in alignments for VEGFA and SPHK2
ENST00000372055strings_subset <- readAAStringSet("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/VEGF_figure/aligns/subset_VEGF_taper.aa")
ENST00000599029strings_subset <- readAAStringSet("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/VEGF_figure/aligns/subset_SPHK2_taper.aa")

##filter out any amino acid not found in any species, which is likely an insertion in a species that was pruned
ENST00000372055string_dat <- tidy_msa(ENST00000372055strings_subset,1,460)
ENST00000372055positions_to_remove <- ENST00000372055string_dat %>%
  group_by(position) %>%
  filter(all(character == "X")) %>%
  distinct(position) %>%
  pull(position)
# Remove those positions
ENST00000372055_filtered <- ENST00000372055string_dat %>%
  filter(!position %in% ENST00000372055positions_to_remove)
# Calculate new total positions
new_total_positions <- nrow(ENST00000372055_filtered) / 11  
# Generate new position sequence
new_positions <- rep(1:new_total_positions, each = 11)
ENST00000372055_filtered$position <- new_positions 
#convert back to AA strings
ENST00000372055_reconstructed <- ENST00000372055_filtered %>%
  group_by(name) %>%                      # Group by species name
  arrange(position) %>%                    # Ensure correct order
  summarise(sequence = paste(character, collapse = ""), .groups = "drop")  # Concatenate characters
# Convert to AAStringSet
AAENST00000372055_reconstructed <- AAStringSet(setNames(ENST00000372055_reconstructed$sequence, ENST00000372055_reconstructed$name))

#and for SPHK2
ENST00000599029string_dat <- tidy_msa(ENST00000599029strings_subset,1,1443)
ENST00000599029positions_to_remove <- ENST00000599029string_dat %>%
  group_by(position) %>%
  filter(all(character == "X")) %>%
  distinct(position) %>%
  pull(position)
# Remove those positions
ENST00000599029_filtered <- ENST00000599029string_dat %>%
  filter(!position %in% ENST00000599029positions_to_remove)
# Calculate new total positions
new_total_positions <- nrow(ENST00000599029_filtered) / 11  
# Generate new position sequence
new_positions <- rep(1:new_total_positions, each = 11)
ENST00000599029_filtered$position <- new_positions
#convert back to AA strings
ENST00000599029_reconstructed <- ENST00000599029_filtered %>%
  group_by(name) %>%                      # Group by species name
  arrange(position) %>%                    # Ensure correct order
  summarise(sequence = paste(character, collapse = ""), .groups = "drop")  # Concatenate characters
# Convert to AAStringSet
AAENST00000599029_reconstructed <- AAStringSet(setNames(ENST00000599029_reconstructed$sequence, ENST00000599029_reconstructed$name))




#boxes and ticks
#########
bandt_frame <- as.data.frame(AAENST00000599029_reconstructed[,1])
bandt_frame$names<-rownames(bandt_frame)
bandt_frame
# Initialize an empty list to store individual data frames
result_list <- list()

# Iterate through each row of the input data frame
for (i in 1:nrow(bandt_frame)) {
  # Extract sequence and species name for the current row
  input_string <- bandt_frame$x[i]
  species_name <- bandt_frame$names[i]
  
  # Convert the string to a vector of characters
  characters <- strsplit(input_string, "")[[1]]
  
  # Create a data frame for the current species
  df <- data.frame(
    position = 1:length(characters),
    character = characters,
    is_X = characters == "X",
    species = species_name
  )
  
  # Initialize variables to store the ranges for the current species
  start_position <- NULL
  end_position <- NULL
  current_value <- NULL
  
  # Initialize an empty data frame to store the result for the current species
  result_df_species <- data.frame(start_position = integer(), end_position = integer(), value = logical(), species = character())
  
  # Iterate through the data frame for the current species
  for (j in 1:nrow(df)) {
    if (is.null(current_value)) {
      # First iteration, initialize variables
      current_value <- df$is_X[j]
      start_position <- df$position[j]
      end_position <- df$position[j]
    } else if (df$is_X[j] == current_value) {
      # Continue the current range
      end_position <- df$position[j]
    } else {
      # End of the current range, add to the result data frame for the current species
      new_row <- data.frame(start_position, end_position + 1, value = current_value, species = species_name)
      result_df_species <- rbind(result_df_species, new_row)
      
      # Reset variables for the next range
      current_value <- df$is_X[j]
      start_position <- df$position[j]
      end_position <- df$position[j]
    }
  }
  
  # Add the last range to the result data frame for the current species
  if (!is.null(current_value)) {
    new_row <- data.frame(start_position, end_position + 1, value = current_value, species = species_name)
    result_df_species <- rbind(result_df_species, new_row)
  }
  
  # Append the result data frame for the current species to the list
  result_list[[i]] <- result_df_species
}
result_list
# Combine individual data frames into a single data frame for all species
result_df <- do.call(rbind, result_list)
result_df
# Print the result data frame
unique(result_df$species)
###
result_df$size <- ifelse(result_df$value == "FALSE", 1, .5)
result_df$end_position1 <- result_df$end_position...1
result_df <- subset(result_df, select = -c(end_position...1))
result_df
#Now this is where things get even more janky
#the bar plots will NOT recognize your species and put them aligned with it
#but will do so based on numerical order in which they appear
#set seq from 1 to how many species you have
#then order of species from bottom branch to top branch
#can be found by continuous printing phylogeny with labeled branch (line 8), ir printing your phylogeny with names but ggtree hates printing words
spec_key <- data.frame(
  A = seq(1, 11, by = 1),
  B = c("sorex_araneus","sus_scrofa","ovis_orientalis","bos_gaurus","macaca_mulatta","homo_sapiens","pan_troglodytes","sciurus_vulgaris","peromyscus_leucopus","mus_musculus","rattus_norvegicus"))
merged_df <- merge(result_df, spec_key, by.x = "species", by.y = "B", all.x = TRUE)
merged_df
order_ENST2 <- as.data.frame(cbind(merged_df$species,merged_df$start_position,merged_df$value, merged_df$size,merged_df$end_position1,merged_df$species,merged_df$A))
order_ENST2$V2<-as.numeric(order_ENST2$V2)
order_ENST2$V3<-as.logical(order_ENST2$V3)
order_ENST2$V4<-as.numeric(order_ENST2$V4)
order_ENST2$V5<-as.numeric(order_ENST2$V5)
order_ENST2$V7<-as.numeric(order_ENST2$V7)




#now need to fix positions on these
VEGFAsnp_data<-read_table("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/VEGF_figure/VEGFA_rfigure_input.txt", 
                          col_names = FALSE)
SPHK2snp_data<-read_table("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/VEGF_figure/SPHK2_rfigure_input.txt", 
                          col_names = FALSE)
#fix by taking positions and subtracting number of positions below them
VEGFAsnp_data <- VEGFAsnp_data %>%
  mutate(
    adjusted_X2 = X2 - sapply(X2, function(pos) sum(ENST00000372055positions_to_remove <= pos))
  )
SPHK2snp_data <- SPHK2snp_data %>%
  mutate(
    adjusted_X2 = X2 - sapply(X2, function(pos) sum(ENST00000599029positions_to_remove <= pos))
  )


## visualize the tree 
p <-ggtree(VEGFjan2025Tree) %<+% VEGFjan2025Tree_images + 
  #geom_tippoint(aes(label=labeling,color=fore),size=3, alpha=.5)+
  geom_tiplab(aes(label=labeling,color=fore), offset=0.04, size =.04)+
  scale_color_manual(values=c("black","lightslateblue"))
p
#note below, False and species change to color of choice
#can have chatGPT whip this up for you quickly by inputing species list and describing what you want, can probably code it with spec_key as well
## and align them based on tree structure
plot_PSG_hyp<- p + geom_facet(panel = "VEGFA Selected Sites", data =order_ENST1, geom = geom_segment,
                              mapping=aes(x=V2,xend=V5,y=V7,yend=V7,size=V4,color = interaction(V3, V6)))+
  geom_facet(panel = "VEGFA Selected Sites", data = VEGFAsnp_data, geom = geom_point, 
             mapping=aes(x = adjusted_X2), shape = '|',color="red", size=5) +
  geom_facet(panel = "SPHK2 Selected Sites", data =order_ENST2, geom = geom_segment,
             mapping=aes(x=V2,xend=V5,y=V7,yend=V7,size=V4,color = interaction(V3, V6))) +
  scale_color_manual(values = species_vector <- c(
    "TRUE.acomys_russatus" = "gray17", "FALSE.acomys_russatus" = "#2F6B8EFF",
    "TRUE.ailuropoda_melanoleuca" = "gray17", "FALSE.ailuropoda_melanoleuca" = "#2F6B8EFF",
    "TRUE.bos_gaurus" = "gray17", "FALSE.bos_gaurus" = "#2F6B8EFF",
    "TRUE.camelus_ferus" = "gray17", "FALSE.camelus_ferus" = "#2F6B8EFF",
    "TRUE.ceratotherium_simum_cottoni" = "gray17", "FALSE.ceratotherium_simum_cottoni" = "#2F6B8EFF",
    "TRUE.cervus_elaphus" = "gray17", "FALSE.cervus_elaphus" = "#2F6B8EFF",
    "TRUE.condylura_cristata" = "gray17", "FALSE.condylura_cristata" = "#2F6B8EFF",
    "TRUE.erinaceus_europaeus" = "gray17", "FALSE.erinaceus_europaeus" = "#2F6B8EFF",
    "TRUE.homo_sapiens" = "gray17", "FALSE.homo_sapiens" = "#2F6B8EFF",
    "TRUE.hylobates_moloch" = "gray17", "FALSE.hylobates_moloch" = "#2F6B8EFF",
    "TRUE.lutra_lutra" = "gray17", "FALSE.lutra_lutra" = "#2F6B8EFF",
    "TRUE.lynx_canadensis" = "gray17", "FALSE.lynx_canadensis" = "#2F6B8EFF",
    "TRUE.macaca_mulatta" = "gray17", "FALSE.macaca_mulatta" = "#2F6B8EFF",
    "TRUE.manis_javanica" = "gray17", "FALSE.manis_javanica" = "#2F6B8EFF",
    "TRUE.martes_zibellina" = "gray17", "FALSE.martes_zibellina" = "#2F6B8EFF",
    "TRUE.microcebus_murinus" = "gray17", "FALSE.microcebus_murinus" = "#2F6B8EFF",
    "TRUE.mus_musculus" = "gray17", "FALSE.mus_musculus" = "#2F6B8EFF",
    "TRUE.mustela_erminea" = "gray17", "FALSE.mustela_erminea" = "lightslateblue",
    "TRUE.mustela_nigripes" = "gray17", "FALSE.mustela_nigripes" = "#2F6B8EFF",
    "TRUE.mustela_putorius_furo" = "gray17", "FALSE.mustela_putorius_furo" = "lightslateblue",
    "TRUE.myotis_myotis" = "gray17", "FALSE.myotis_myotis" = "#2F6B8EFF",
    "TRUE.oryctolagus_cuniculus" = "gray17", "FALSE.oryctolagus_cuniculus" = "#2F6B8EFF",
    "TRUE.ovis_orientalis" = "gray17", "FALSE.ovis_orientalis" = "#2F6B8EFF",
    "TRUE.pan_troglodytes" = "gray17", "FALSE.pan_troglodytes" = "#2F6B8EFF",
    "TRUE.papio_anubis" = "gray17", "FALSE.papio_anubis" = "#2F6B8EFF",
    "TRUE.peromyscus_leucopus" = "gray17", "FALSE.peromyscus_leucopus" = "#2F6B8EFF",
    "TRUE.phocoena_sinus" = "gray17", "FALSE.phocoena_sinus" = "#2F6B8EFF",
    "TRUE.phyllostomus_discolor" = "gray17", "FALSE.phyllostomus_discolor" = "#2F6B8EFF",
    "TRUE.pipistrellus_pipistrellus" = "gray17", "FALSE.pipistrellus_pipistrellus" = "#2F6B8EFF",
    "TRUE.rattus_norvegicus" = "gray17", "FALSE.rattus_norvegicus" = "#2F6B8EFF",
    "TRUE.rhinolophus_ferrumequinum" = "gray17", "FALSE.rhinolophus_ferrumequinum" = "#2F6B8EFF",
    "TRUE.scalopus_aquaticus" = "gray17", "FALSE.scalopus_aquaticus" = "#2F6B8EFF",
    "TRUE.sciurus_vulgaris" = "gray17", "FALSE.sciurus_vulgaris" = "#2F6B8EFF",
    "TRUE.solenodon_paradoxus" = "gray17", "FALSE.solenodon_paradoxus" = "#2F6B8EFF",
    "TRUE.sorex_araneus" = "gray17", "FALSE.sorex_araneus" = "lightslateblue",
    "TRUE.suncus_etruscus" = "gray17", "FALSE.suncus_etruscus" = "lightslateblue",
    "TRUE.sus_scrofa" = "gray17", "FALSE.sus_scrofa" = "#2F6B8EFF",
    "TRUE.talpa_occidentalis" = "gray17", "FALSE.talpa_occidentalis" = "#2F6B8EFF",
    "TRUE.tupaia_chinensis" = "gray17", "FALSE.tupaia_chinensis" = "#2F6B8EFF",
    "TRUE.vulpes_lagopus" = "gray17", "FALSE.vulpes_lagopus" = "#2F6B8EFF"
  ))+
  geom_facet(panel = "SPHK2 Selected Sites", data = SPHK2snp_data, geom = geom_point, 
             mapping=aes(x = adjusted_X2), shape = '|',color="red", size=5) +
  theme_tree2(legend.position = 'none')
facet_widths(plot_PSG_hyp,c(.2,.42,1.407))
p
gheatmap(p, hyp_heat2,offset = 2,width = 6,colnames_angle=270,hjust=0,colnames_offset_y =.3)+
  scale_fill_viridis_c(option="A", name="continuous\nvalue")
