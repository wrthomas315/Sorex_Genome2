###Graphics for genes
ENST00000378126_tree <- read.tree("~//Sorex_Genome2/analysis/5_expression/1_for_visualizations/PCDHA6/noFSENST00000378126.macse.aligned.fa.taper.fa.nh")
#phylopic depreciated but still run to create needed frame for fore
ENST00000378126_tree_images <- ggimage::phylopic_uid(ENST00000378126_tree$tip.label)
ENST00000378126_tree_images$labeling <-  ENST00000378126_tree_images$name
ENST00000378126_tree_images$labeling
#c("macaca_mulatta","hylobates_moloch","homo_sapiens","pan_troglodytes","sciurus_vulgaris","peromyscus_leucopus","acomys_russatus","mus_musculus","rattus_norvegicus","talpa_occidentalis","suncus_etruscus","sorex_araneus","rhinolophus_ferrumequinum","phyllostomus_discolor","pipistrellus_pipistrellus","myotis_myotis","manis_javanica","lutra_lutra","mustela_erminea","mustela_putorius_furo","camelus_ferus","sus_scrofa","phocoena_sinus","ovis_orientalis","bos_gaurus")
ENST00000378126_tree_images$fore <- c(rep("Background",1),rep("Foreground",2),rep("Background",10),rep("Foreground",2),rep("Background",10))
#ENST00000378126_tree_images$fore <- c(rep("Background",22),rep("Foreground",1),rep("Background",2))
#check to make sure foreground and background are w/ correct species
ENST00000378126_tree_images
ENST00000378126_phy <-ggtree(ENST00000378126_tree,branch.length="none") %<+% ENST00000378126_tree_images + 
  geom_tippoint(aes(label=labeling,color=fore),size=3, alpha=.5)+
  geom_tiplab(aes(label=labeling,color=fore), offset=0.04, size =.04)+
  scale_color_manual(values=c("black","lightslateblue"))
ENST00000378126_phy
ENST00000378126strings <- readAAStringSet("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/PCDHA6/ENST00000378126.aa")
#now where do we want to look at in this gene
ENST00000378126string_dat <- tidy_msa(ENST00000378126strings,105,140)
msa_ENST00000378126 <-ENST00000378126_phy + geom_facet(geom = geom_msa, data = ENST00000378126string_dat,  panel = 'msa')+theme_tree2(legend.position = 'none')
facet_widths(msa_ENST00000378126,c(.10,.85))

#boxes and ticks
#########
bandt_frame <- as.data.frame(ENST00000378126strings[,1])
bandt_frame$names<-rownames(bandt_frame)
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

# Combine individual data frames into a single data frame for all species
result_df <- do.call(rbind, result_list)

# Print the result data frame
print(result_df)
unique(result_df$species)
###
result_df$size <- ifelse(result_df$value == "FALSE", 1, .5)
result_df$end_position1 <- result_df$end_position...1
result_df <- subset(result_df, select = -c(end_position...1))
spec_key <- data.frame(
  A = seq(1, 25, by = 1),
  B = c("macaca_mulatta","hylobates_moloch","homo_sapiens","pan_troglodytes","sciurus_vulgaris","peromyscus_leucopus","acomys_russatus","mus_musculus","rattus_norvegicus","talpa_occidentalis","suncus_etruscus","sorex_araneus","rhinolophus_ferrumequinum","phyllostomus_discolor","pipistrellus_pipistrellus","myotis_myotis","manis_javanica","lutra_lutra","mustela_erminea","mustela_putorius_furo","camelus_ferus","sus_scrofa","phocoena_sinus","ovis_orientalis","bos_gaurus"))
result_df
merged_df <- merge(result_df, spec_key, by.x = "species", by.y = "B", all.x = TRUE)
merged_df
order_ENST2 <- as.data.frame(cbind(merged_df$species,merged_df$start_position,merged_df$value, merged_df$size,merged_df$end_position1,merged_df$species,merged_df$A))
order_ENST2$V2<-as.numeric(order_ENST2$V2)
order_ENST2$V3<-as.logical(order_ENST2$V3)
order_ENST2$V4<-as.numeric(order_ENST2$V4)
order_ENST2$V5<-as.numeric(order_ENST2$V5)
order_ENST2$V7<-as.numeric(order_ENST2$V7)
unique(order_ENST2$V1)
# Plot using ggplot2 with just ticks
str(result_df)
str(order_ENST)
ggplot(data = order_ENST2) +
  geom_segment(aes(x = V1, xend = V1, y = V2, yend = V5, color = interaction(V3, V1), size = V4)) +
  coord_flip() +
  ylab("Value") +
  scale_color_manual(values = c("FALSE.homo_sapiens" = "pink", "TRUE.homo_sapiens" = "black", "FALSE.macaca_mulatta" = "black", "TRUE.macaca_mulatta" = "black", "FALSE.pan_troglodytes" = "blue", "TRUE.pan_troglodytes" = "black")) +
  scale_size_identity() +
  theme_minimal()
PCDHA6snp_data<-read_table("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/PCDHA6PCDHA6_rfigure_input.txt", 
                           col_names = FALSE)
## visualize the tree 
p <-ggtree(ENST00000378126_tree) %<+% ENST00000378126_tree_images + 
  #geom_tippoint(aes(label=labeling,color=fore),size=3, alpha=.5)+
  geom_tiplab(aes(label=labeling,color=fore), offset=0.04, size =.04)+
  scale_color_manual(values=c("black","lightslateblue"))

## and align them based on tree structure
plot_PCDHA6<-p + geom_facet(panel = "PCDHA6 Selected Sites", data =order_ENST2, geom = geom_segment,
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
  geom_facet(panel = "PCDHA6 Selected Sites", data = PCDHA6snp_data, geom = geom_point, 
             mapping=aes(x = X2), shape = '|',color="red", size=5) +
  theme_tree2(legend.position = 'none')
#########
facet_widths(plot_PCDHA6,c(.1,1.6))
