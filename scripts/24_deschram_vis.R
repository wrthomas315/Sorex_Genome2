########
###SYNTENY FIGURES
#load library
library(syntenyPlotteR)
###Make figures for pairwise comparison to APCF
#Some versions of reformatting and drawing commented out, but left in if you want figure references and reversal to change
#reformat data
#reformat.syntenyData("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/APCF_sorAra.merged.map.mod","~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/APCF_sorAra.merged.map.reformat", reference.species = "APCF", target.species = "sorAra")
reformat.syntenyData("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/APCF_sorAra.merged.map.mod","~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/reformat/APCF_sorAra.merged.map.reformat",reference.species = "APCF", target.species = "sorAra")
reformat.syntenyData("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/APCF_sorAra.merged.map.mod","~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/reformat/revAPCF_sorAra.merged.map.reformat",reference.species = "sorAra", target.species = "APCF")
#reformat swithcing sor to reference
#draw linear from both aspects
draw.linear("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/figures/pairwise", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/reformat/APCF_size.txt", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/reformat/APCF_sorAra.merged.map.reformat.txt")
draw.linear("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/figures/revpairwise", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/reformat/revAPCF_size.txt", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/reformat/revAPCF_sorAra.merged.map.reformat.txt")
#and ideogram
draw.ideogram("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/APCF_sorAra.merged.map.reformat.txt", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/rev.sorAPCF.combo.chromsize", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/figures/ideopairwise")
draw.ideogram("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/revAPCF_sorAra.merged.map.reformat.txt", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/sorAPCF.combo.chromsize", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/APCF/figures/revideopairwise",h=13)

###Make figures for pairwise comparison to APCF
reformat.syntenyData("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/SFs/SFs_sunEtr/Orthology.Blocks.mod22","~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/SFs/SFs_sunEtr/Orthology.Blocks_sun.reformat", reference.species = "sorAra", target.species = "sunEtr")
#reformat with sor as reference
#draw linear from both aspects
draw.linear("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/figures/sunpairwise", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunEtr/OFSsorsun.combo.chromsize", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/SFs/SFs_sunEtr/Orthology.Blocks_sun.reformat.txt")
#draw.linear("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunEtr/figures/revpairwise", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunEtr/revOFSsorsun.combo.chromsize", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunEtr/revOrthology.Blocks_sun.reformat.txt")
#maybe the ideogram will look better
draw.ideogram("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/SFs/SFs_sunEtr/Orthology.Blocks_sun.reformat.txt", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunEtr/OFSsorsun.combo.chromsize", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/shrew_wEur/shrew_Eul_1000/figures/sunideopairwise")


###figure combining both sor sun and APCF
reformat.syntenyData("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunAPCF/APCF_sunEtr.merged.map22","~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunAPCF/APCF_sunEtr.merged.map22.reformat", reference.species = "APCF", target.species = "sunEtr")
draw.linear("~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunAPCF/figures/pairwise", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunAPCF/rev.sunAPCFsor.combo.chromsize", "~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunEtr/Orthology.Blocks_sun.reformat.txt","~/Sorex_Genome2/analysis/6_deschrambler/1_shrew_Eul_1000/SV/PU/sunAPCF/APCF_sunEtr.merged.map22.reformat.txt")