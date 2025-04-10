#upset plot from saved file if ensembl not working (it tends to not to)
library(readr)
library(dplyr)
library(tidyr)
hyphytotal <- read_delim("~/Sorex_Genome2/analysis/3_hyphy/hyphytotal.txt", 
                         delim = "\t", escape_double = FALSE)

PSGsorAra<-subset(hyphytotal,hyphytotal$X3 == "sorex_araneus")
PSGmusMus<-subset(hyphytotal,hyphytotal$X3 == "mus_musculus")
PSGmusErm<-subset(hyphytotal,hyphytotal$X3 == "mustela_erminea")
PSGsolPar<-subset(hyphytotal,hyphytotal$X3 == "solenodon_paradoxus")
PSGscaAqu<-subset(hyphytotal,hyphytotal$X3 == "scalopus_aquaticus")
PSGeriEur<-subset(hyphytotal,hyphytotal$X3 == "erinaceus_europaeus")
PSGconCri<-subset(hyphytotal,hyphytotal$X3 == "condylura_cristata")
PSGmusPutFur<-subset(hyphytotal,hyphytotal$X3 == "mustela_putorius_furo")
PSGmarZib<-subset(hyphytotal,hyphytotal$X3 == "martes_zibellina")
PSGlutLut<-subset(hyphytotal,hyphytotal$X3 == "lutra_lutra")
PSGmusNig<-subset(hyphytotal,hyphytotal$X3 == "mustela_nigripes")
PSGailMel<-subset(hyphytotal,hyphytotal$X3 == "ailuropoda_melanoleuca")
PSGlynCan<-subset(hyphytotal,hyphytotal$X3 == "lynx_canadensis")
PSGvulVul<-subset(hyphytotal,hyphytotal$X3 == "vulpes_lagopus")
PSGmanJav<-subset(hyphytotal,hyphytotal$X3 == "manis_javanica")
PSGcerSim<-subset(hyphytotal,hyphytotal$X3 == "ceratotherium_simum_cottoni")
PSGsusScr<-subset(hyphytotal,hyphytotal$X3 == "sus_scrofa")
PSGphoSin<-subset(hyphytotal,hyphytotal$X3 == "phocoena_sinus")
PSGoviOri<-subset(hyphytotal,hyphytotal$X3 == "ovis_orientalis")
PSGbosGau<-subset(hyphytotal,hyphytotal$X3 == "bos_gaurus")
PSGcerEla<-subset(hyphytotal,hyphytotal$X3 == "cervus_elaphus")
PSGcamFer<-subset(hyphytotal,hyphytotal$X3 == "camelus_ferus")
PSGrhiFer<-subset(hyphytotal,hyphytotal$X3 == "rhinolophus_ferrumequinum")
PSGmyoMyo<-subset(hyphytotal,hyphytotal$X3 == "myotis_myotis")
PSGphyDis<-subset(hyphytotal,hyphytotal$X3 == "phyllostomus_discolor")
PSGpipPip<-subset(hyphytotal,hyphytotal$X3 == "pipistrellus_pipistrellus")
PSGoryCun<-subset(hyphytotal,hyphytotal$X3 == "oryctolagus_cuniculus")
PSGsciVul<-subset(hyphytotal,hyphytotal$X3 == "sciurus_vulgaris")
PSGratNor<-subset(hyphytotal,hyphytotal$X3 == "rattus_norvegicus")
PSGperLeu<-subset(hyphytotal,hyphytotal$X3 == "peromyscus_leucopus")
PSGacoRus<-subset(hyphytotal,hyphytotal$X3 == "acomys_russatus")
PSGmacMul<-subset(hyphytotal,hyphytotal$X3 == "macaca_mulatta")
PSGpanTro<-subset(hyphytotal,hyphytotal$X3 == "pan_troglodytes")
PSGhylMol<-subset(hyphytotal,hyphytotal$X3 == "hylobates_moloch")
PSGpapAnu<-subset(hyphytotal,hyphytotal$X3 == "papio_anubis")
PSGmicMur<-subset(hyphytotal,hyphytotal$X3 == "microcebus_murinus")
PSGtupChi<-subset(hyphytotal,hyphytotal$X3 == "tupaia_chinensis")
PSGtalOcc<-subset(hyphytotal,hyphytotal$X3 == "talpa_occidentalis")
PSGsunEtr<-subset(hyphytotal,hyphytotal$X3 == "suncus_etruscus")
PSGhomSap<-subset(hyphytotal,hyphytotal$X3 == "homo_sapiens")

combined_vector <- rbind(PSGhomSap,PSGtalOcc,PSGtupChi,PSGmicMur,PSGpapAnu,PSGhylMol,PSGpanTro,PSGmacMul,PSGacoRus,PSGperLeu,PSGratNor,
                         PSGsciVul,PSGoryCun,PSGpipPip,PSGphyDis,PSGmyoMyo,PSGrhiFer,PSGcamFer,PSGcerEla,PSGbosGau,PSGoviOri,PSGphoSin,
                         PSGsusScr,PSGcerSim,PSGmanJav,PSGvulVul,PSGlynCan,PSGailMel,PSGlutLut,PSGmarZib,PSGconCri,PSGeriEur,PSGscaAqu,
                         PSGsolPar,PSGmusMus,PSGmusNig)

# Iterate over rows in absrel results
GMenst_counts <- table(combined_vector$X2)
GMresult_gt5 <- names(GMenst_counts[GMenst_counts > 4])
GNlooseCONV10 <-unique(c(setdiff(Reduce(intersect, list(PSGsunEtr$X2, PSGsorAra$X2,PSGmusErm$X2)),GMresult_gt5),setdiff(Reduce(intersect, list(PSGmusPutFur$X2, PSGsorAra$X2,PSGmusErm$X2)),GMresult_gt5),setdiff(Reduce(intersect, list(PSGsunEtr$X2, PSGsorAra$X2,PSGmusPutFur$X2)),GMresult_gt5)))
#1 (PCDHA6) convergent strict background
setdiff(Reduce(intersect, list(PSGmusErm$X2, PSGsorAra$X2,PSGmusPutFur$X2)),names(GMenst_counts[GMenst_counts > 0]))
#table 2 check
unique(c(setdiff(Reduce(intersect, list(PSGsunEtr$X2, PSGsorAra$X2,PSGmusPutFur$X2)),GMresult_gt5)))
#
#Final upset for supplement
unique(combined_vector$X2)
GMresult_gt5 <- names(GMenst_counts[GMenst_counts > 4])
gene_listUPSET8 <- list(Mustela_erminea = PSGmusErm$X2,Sorex_araneus = PSGsorAra$X2,Mustela_putorius_furo = PSGmusPutFur$X2,Suncus_etruscus = PSGsunEtr$X2,Loose =names(GMenst_counts[GMenst_counts < 5]), Strict=unique(combined_vector$X2))
upset_gene_listUPSET8 <- fromList(gene_listUPSET8)
upset(upset_gene_listUPSET8,nsets = 6,order.by = "freq",)



#Permutation test
incDehn_combined_vector <-rbind(combined_vector, PSGsorAra, PSGmusErm, PSGmusPutFur, PSGsunEtr)
totgeneList <- read_csv("~/Sorex_Genome2/analysis/3_hyphy/geneList.txt", 
                        col_names = FALSE)
totgeneList_format<- as.data.frame(totgeneList$X1)
colnames(totgeneList_format) <- c("GeneID")
incDehn_combined_format <- as.data.frame(cbind(incDehn_combined_vector$X2, incDehn_combined_vector$X3))
colnames(incDehn_combined_format) <- c("GeneID", "Species")
merged_data <- totgeneList_format %>%
  left_join(incDehn_combined_format, by = "GeneID") %>%
  mutate(Species = ifelse(is.na(Species), "none", Species))
merged_data
binary_matrix <- merged_data %>%
  mutate(present = ifelse(Species == "none", 0, 1)) %>%
  select(GeneID, Species, present) %>%
  pivot_wider(names_from = Species, values_from = present, values_fill = 0)

# Step 3: Rename 'none' column to represent species that did not show selection
gene_data <- binary_matrix %>% select(-none)

# View the binary matrix
print(gene_data)
#compare look good
colSums(gene_data[ , -1])


# Define the columns you're interested in for selection counts
selection_columns <- c("mustela_putorius_furo", "suncus_etruscus", "mustela_erminea")
remaining_columns <- setdiff(names(gene_data), c("GeneID", "sorex_araneus", selection_columns))

# Apply the filtering conditions
filtered_genes <- gene_data %>%
  filter(
    sorex_araneus == 1,  # Condition 1: Positive selection in sorex_araneus
    rowSums(select(., selection_columns) == 1) >= 2,  # Condition 2: Positive selection in at least 2 of the 3 specified columns
    rowSums(select(., remaining_columns) == 1) <= 4  # Condition 3: Not under positive selection in > 4 of the remaining columns (excluding sorex_araneus)
  )

# Count the number of qualifying genes
count_filtered_genes <- nrow(filtered_genes)
count_filtered_genes

#permute
calc_qualifying_genes <- function(df) {
  filtered_genes <- df %>%
    filter(
      sorex_araneus == 1,  # Condition 1: Positive selection in sorex_araneus
      rowSums(select(., selection_columns) == 1) >= 2,  # Condition 2: Positive selection in at least 2 of the 3 specified columns
      rowSums(select(., remaining_columns) == 1) <= 4  # Condition 3: Not under positive selection in > 4 of the remaining columns
    )
  return(nrow(filtered_genes))
}
#shuffle function
shuffle_species_column <- function(df, species_column) {
  # Shuffle the selection values (1's and 0's) for a single species
  shuffled_values <- sample(df[[species_column]])
  df[[species_column]] <- shuffled_values
  return(df)
}
sum(sample(gene_data[["mustela_erminea"]]))
#permute function
permute_and_test <- function(gene_data, num_permutations) {
  permuted_results <- numeric(num_permutations)
  
  for (i in 1:num_permutations) {
    permuted_df <- gene_data
    
    # Shuffle the species columns
    for (species in selection_columns) {
      permuted_df <- shuffle_species_column(permuted_df, species)
    }
    
    # Calculate the number of qualifying genes for this permutation
    permuted_results[i] <- calc_qualifying_genes(permuted_df)
  }
  
  return(permuted_results)
}
permuted_results_res <- permute_and_test(gene_data, 100000)
permuted_results_res
#pvalue
observed_count <- count_filtered_genes
p_value <- mean(permuted_results_res >= observed_count)
p_value


# Step 7: Plotting the permutation results
hist(permuted_results_res, breaks = 50, main = "Permutation Test Results",
     xlab = "Number of Convergent Genes", col = "lightblue", border = "black",
     xlim =c(0,8), xaxs = "i", yaxs = "i")
abline(v = observed_convergence, col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Observed", ""), col = "red", lwd = 1.5, lty = 2)
