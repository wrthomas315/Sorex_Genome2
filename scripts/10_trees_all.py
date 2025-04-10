from ete3 import Tree
import re
import os

species_list = []
#species_list.append = []

for filename in os.listdir("~/Sorex_Genome2/data/8_taper_align/"):
	if filename.endswith(".fa"):	
		file= open("~/Sorex_Genome2/data/8_taper_align/" + filename,"r")
		string = file.read()
		species = string.split()
		del species[1::2]
		new = [s.replace('>', '') for s in species]
		file.close()
		treefile = open("~/Sorex_Genome2/data/9_trees/Alvarez_Carretero2022_4075sp.pruned.nh","r")
		tree = treefile.read()
		t = Tree(tree, format = 1)
		t.prune(new)
		t.unroot()
		t.write(format = 1, outfile = "~/Sorex_Genome2/data/9_trees/unrooted/" + filename + ".nh")
