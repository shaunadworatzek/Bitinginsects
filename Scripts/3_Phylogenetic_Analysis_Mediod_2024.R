#Data required is KBIMP2024_filteredCOI.tsv from the data cleaning step and 
#the outgroup file 
#the output files will be the updated species file and the two 
#phylogenies which are PNGs 

library(stringr)
library(tidyverse)
library(readr)
library(viridis)
library(Biostrings) 
library(ape)     
library(muscle) 
library(phangorn) 
library(parallel) 
library(ggplot2)    
library(ggtree)     
library(seqinr)
library(DECIPHER)


KBIMP2024 <- read_tsv(file = "processed-data/KBIMP2024_filteredCOI.tsv")
Outgroup <- read_csv(file = "raw-data2/Outgroup.csv")
#### PART 1 - Phylogenetic analysis and resolution of undefined Mosquito species ---- 

#filtering for Culicidae, and defining the root as a chiromidae sequence in the data set 

KBIMP2024 <- KBIMP2024 %>%
  bind_rows(Outgroup)

kbimp_mos_DNA_df <- KBIMP2024 %>%
  filter(Family %in% c("Culicidae", "Outgroup")) %>%
  select(Sample, Sequence) %>%
  
  

#preforming the alignment and converting to a PhyDat

kbimp_mos_DNA <- DNAStringSet(kbimp_mos_DNA_df$Sequence)

names(kbimp_mos_DNA) <- kbimp_mos_DNA_df$Sample

alighned_kbimpmos_DNA <- DNAStringSet(muscle::muscle(kbimp_mos_DNA))

kbimp_mos_phydat <- as.phyDat(alighned_kbimpmos_DNA, type = "DNA")
class(kbimp_mos_phydat) # is a "phyDat" object
length(kbimp_mos_phydat) # 201 has the  species as seen before 

#here before I do model testing I would do a visualization step to visualize the alignment 

###### determining medoid sequences ######

modelTest <- modelTest(kbimp_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 
#JC is the best model because it has the lowest AIC, second lowest BIC and high loglik

dist.kbimp.mos <- dist.ml(kbimp_mos_phydat, ratio = TRUE, model = "JC69") 

#cluster sequences

hc <- hclust(dist.kbimp.mos)

groups <- cutree(hc, k = 15) 
#done multiple times to determine the optimal number of groups based on my data - 
#What criteria did I use to determine the optimum number of clusters 
#- internal criteria of cluster cohesion - external measures like by comparing to species labels 
#Do the internal test to show this because it can be the case that some species have more 
#variability and the 


table(groups) #number of sequences in each group

cluster_sizes <- table(groups)

medoids_2024 <- sapply(unique(groups), function(g) {
  members <- names(groups[groups == g])
  find_medoid(dist.kbimp.mos, members)
})
#Identifying a single sequence which is the most central of the group instead 
#of making a new average sequence - If I were defining species grouping de novo - there is an argument to make to be that the 
#equivalent of blasting to each individual to each individual unique sequence - an alternative approach would be filtering to unique sequences 
#and then comparing unique sequences to a data base - when I am assigning biological data I am relying on using species names and studies of vector status - don't delete just compare and see if I get consistent results across the two at least in appendex how many total sequences do I have, how many total sequences do I have - what the min and max divergence - not expecting 10-15% within a species probably 2-4% - cite numbers from the literature. 

medoid_phy_2024 <- kbimp_mos_phydat[medoids_2024]

medoid_phy_2024 

###### Building the tree ######

dist.medoid.mos <- dist.ml(medoid_phy_2024, ratio = TRUE, model = "JC69") 

#neighbor joining method

NJtree.kbimp.mos <- NJ(dist.medoid.mos)

plot(NJtree.kbimp.mos)

length(NJtree.kbimp.mos$tip.label) #1272 tips as expected 

# Fit the initial tree using a simple pml 

pml.tree.kbimp.mos <- pml(NJtree.kbimp.mos, kbimp_mos_phydat, k = 4, model = "GTR+I", method = "unrooted")

plot(pml.tree.kbimp.mos$tree) 

pml.tree.kbimp.mosop <- pml_bb(pml.tree.kbimp.mos, model = "GTR+I")

#bootstrapping to determine reliability

bs <- bootstrap.pml(pml.tree.kbimp.mosop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bs <- plotBS(pml.tree.kbimp.mos$tree, bs)

#rooting the tree

rooted.bstree.mos <- root(tree_with_bs, outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.mos)

tree_with_bs <- plotBS(rooted.bstree.mos, bs)

###### Species assignment of unknown medoid groups ######

#After BLASTing these sequences I change the name of these 14 sequences and add the species name 

speciesdata <- KBIMP2024%>%
  filter(Sample %in% c("Outgroup",
    "KBIMP-_005_A10" ,"KBIMP-_005_B1", 
     "KBIMP-_005_A2" ,"KBIMP-_005_C7" ,
   "KBIMP-_005_D8" , "KBIMP-_006_B4", 
   "KBIMP-_006_G1",  "KBIMP_004_F7" , 
   "KBIMP-_008_A1",  "KBIMP-_007_F10",
    "KBIMP-_007_F2" , "KBIMP-_007_H11",
    "KBIMP-_007_H4" , "KBIMP-_008_E3" ,
    "KBIMP-_008_F3" )) %>%
  
  select(Sample, Species) %>%

  mutate(ID = NA_character_) %>%
  
  mutate(Species = if_else(Sample == "KBIMP-_005_B1", "Aedes nigripes/impiger", Species)) %>% #92.31 nigripes on bold, 100% on NCBI
  mutate(Species = if_else(Sample == "KBIMP_004_F7", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_006_G1", "Aedes nigripes/impiger", Species)) %>% #100% impeger on NCBI slightly less sure on BOLD
  mutate(Species = if_else(Sample == "KBIMP-_005_D8", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_C6", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_006_B4", "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_A1", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_F2", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_E3", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_H4", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_H11", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_003_G5", "Aedes punctor/Aedes hexodontus", Species)) %>%
  
  mutate(ID = if_else(Sample %in% c(
    "KBIMP-_008_A1",
    "KBIMP-_005_B1",
    "KBIMP_004_F7",
    "KBIMP-_006_G1",
    "KBIMP-_005_D8",
    "KBIMP-_007_C6",
    "KBIMP-_006_B4",
    "KBIMP-_007_G10",
    "KBIMP-_007_F2",
    "KBIMP-_008_E3",
    "KBIMP-_007_H4",
    "KBIMP-_007_H11",
    "KBIMP_003_G5"
  ), "Sequence similarity", ID)) %>%

  mutate(ID = if_else(is.na(ID), "Probabilistic taxonomic assignment", ID))


##### Finalized Figure ######

# match trait data and tree data

tree_tips_mos <- rooted.bstree.mos$tip.label 

sample_name_mos <- speciesdata$Sample 

# checking if all species names are present in the tree

all(tree_tips_mos %in% sample_name_mos) 

all(sample_name_mos %in% tree_tips_mos)

setdiff(tree_tips_mos, sample_name_mos)

setdiff(sample_name_mos, tree_tips_mos)

speciesdata <- speciesdata[match(tree_tips_mos, speciesdata$Sample), ] %>%
  rownames_to_column(var = "cluster")

speciesdata$cluster_size <- cluster_sizes[as.character(speciesdata$cluster)]

group_colors <-  c("Sequence similarity" = "#5050FFFF", "Probabilistic taxonomic assignment" = "#CE3D32FF")

speciesdata$color <- group_colors[speciesdata$ID]

#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.mos)

speciesdata$tip_label_new <- paste0(speciesdata$Species," (", speciesdata$cluster_size, ")")

#produce the tree

mosquitotree2024 <- (
  ggtree(rooted.bstree.mos, layout = "rectangular", branch.length = TRUE) +
    
    geom_tippoint(aes(color = speciesdata$ID[node]), size = 3) +
    
    scale_color_manual(values = group_colors, name = "Identification Method") +
    
    geom_text(aes(label = speciesdata$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
    
    theme(legend.title = element_text(size = 14),      
          legend.text = element_text(size = 12), legend.position = "top")
)

plot(mosquitotree2024)

#save the tree

ggsave("plots/mosphylomoregroups.png", plot = mosquitotree2024, width = 10, height = 4, dpi = 300)




seqs_punc1_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 1]]
seqs_nig_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 2]]
seqs_punc2_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 3]]
seqs_punc3_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 4]]
seqs_imp2_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 5]]
seqs_imp3_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 6]]
seqs_punc4_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 7]]
seqs_ex_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 8]]
#9 is alaskenesis 
seqs_ex2_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 10]]
seqs_punc5_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 11]]
seqs_punc6_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 12]]
seqs_ex3_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 13]]
#14 is inornata

KBIMP_updatedspecies_mos_2024 <- KBIMP2024 %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punc1_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punc2_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punc3_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punc4_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punc5_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punc6_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_nig_mos) , "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_imp2_mos) , "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_imp3_mos) , "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_ex_mos) , "Aedes excrucians", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_ex2_mos) , "Aedes excrucians", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_ex3_mos) , "Aedes excrucians", Species)) 

#### - PART 2 - Black flies - Phylogenetic analysis and resolution of undefined species ----

#filtering for Simuliidae, and defining the root as a chirmoidae sequence in the data set 

kbimp_bf_DNA_df <- KBIMP2024%>%
  mutate(Family = if_else(Sample == "KBIMP-_005_A2", "Root", Family)) %>%
  filter(Family %in% c("Simuliidae", "Outgroup")) %>%
  mutate(Sample = if_else(Sample == "KBIMP-_005_A2", "Root", Sample)) %>%
  select(Sample, Sequence)  %>%
  bind_rows(Outgroup)

# Convert the 'sequence' column to a DNAStringSet

kbimp_bf_DNA <- DNAStringSet(kbimp_bf_DNA_df$Sequence)

names(kbimp_bf_DNA) <- kbimp_bf_DNA_df$Sample

#preforming an alignment 

alighned_kbimpbf_DNA <- DNAStringSet(muscle::muscle(kbimp_bf_DNA))

#converting fasta file into a phyDat format 

kbimp_bf_phydat <- as.phyDat(alighned_kbimpbf_DNA, type = "DNA")

class(kbimp_bf_phydat) 

length(kbimp_bf_phydat) 

#determining the best model for the ML tree 

#modelTest <- modelTest(KGLGTK_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) # JC is the best model because it has the lowest AIC, second lowest BIC and high loglik

# creating a distance matrix

dist.kbimp.bf <- dist.ml(kbimp_bf_phydat, ratio = TRUE, model = "JC69") 

# Hierarchical clustering

hc_bf <- hclust(dist.kbimp.bf)

# get cluster memberships

groups <- cutree(hc_bf, k = 20) 
# I did this multiple different times to determine the optimal number of groups based on my data 

table(groups) 
cluster_sizes <- table(groups) 

#find mediods 

medoids <- sapply(unique(groups), function(g) {
  members <- names(groups[groups == g])
  find_medoid(dist.kbimp.bf, members)})

medoids

# convert subset to a phyDat object

medoid_phy <- kbimp_bf_phydat[medoids]

medoid_phy

#create a new dist matrix

dist.medoid.bf <- dist.ml(medoid_phy, ratio = TRUE, model = "JC69") 

# creating a tree using the neighbor joining method

NJtree.kbimp.bf <- NJ(dist.medoid.bf)

plot(NJtree.kbimp.bf)

length(NJtree.kbimp.bf$tip.label) 

# Fit the initial tree using a simple pml

pml.tree.kbimp.bf <- pml(NJtree.kbimp.bf, kbimp_bf_phydat, k = 4, model = "GTR+I", method = "unrooted")

plot(pml.tree.kbimp.bf$tree) 

pml.tree.kbimp.bfop <- pml_bb(pml.tree.kbimp.bf, model = "GTR+I")

#bootstrapping analysis for tree

bs.bf <- bootstrap.pml(pml.tree.kbimp.bfop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bsbf <- plotBS(pml.tree.kbimp.bf$tree, bs.bf)

#rooting the 

rooted.bstree.bf <- root(tree_with_bsbf, outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.bf)

tree_with_bs.bf <- plotBS(rooted.bstree.bf, bs.bf)

#After BLASTing these sequences I change the name of these 14 sequences and add the species name 

speciesdata_bf <- KBIMP2024%>%
  
  filter(Sample %in% c("Outgroup", "KBIMP-_005_G1", "KBIMP-_007_B9",
                       "KBIMP-_006_C5", "KBIMP_003_F10", "KBIMP-_005_A1", 
                       "KBIMP-_008_D5", "KBIMP-_008_D4", "KBIMP-_005_A12",
                       "KBIMP-_008_B2", "KBIMP-_005_F4" , "KBIMP_003_E7", 
                       "KBIMP-_007_A1", "KBIMP-_007_E10", "KBIMP_003_C5" ,
                       "KBIMP_003_F10" , "KBIMP_003_B12", "KBIMP-_007_C5", 
                       "KBIMP-_007_C11", "KBIMP-_008_A5", "KBIMP_003_F2")) %>%
  
  select(Sample, Species) %>%
  
  mutate(ID = NA_character_) %>%
  
  mutate(Species = if_else(Sample == "KBIMP-_007_B9", "Simulium congareenarum", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_003_F10", "Metacnephia borealis", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_D4", "Simulium craigi", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_005_A12", "Simulium sp.", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_B2", "Simulium subpusillum", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_E10", "Simulium tuberosum", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_C11", "Simulium arcticum complex", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_A5", "Simulium venustum complex", Species)) %>%
  
  mutate(ID = if_else(Sample %in% c(
    "KBIMP-_007_B9", "KBIMP_003_F10", "KBIMP-_008_D4",
    "KBIMP-_005_A12", "KBIMP-_007_E10", 
    "KBIMP-_007_C11",  "KBIMP-_008_A5",
    "KBIMP-_008_B2"), "Sequence similarity", ID)) %>%
  
  mutate(ID = if_else(is.na(ID), "Probabilistic taxonomic assignment", ID))


#making a tree with the sampling locations labels on the tree tips 

tree_tips_bf <- rooted.bstree.bf$tip.label

sample_name_bf <- speciesdata_bf$Sample 

all(tree_tips_bf %in% sample_name_bf) 

all(sample_name_bf %in% tree_tips_bf) 

setdiff(tree_tips_bf, sample_name_bf)

setdiff(sample_name_bf, tree_tips_bf)

#preping tip metadata for plotting 

speciesdata_bf <- speciesdata_bf[match(tree_tips_bf, speciesdata_bf$Sample), ] %>%
  rownames_to_column(var = "cluster")

speciesdata_bf$cluster_size <- cluster_sizes[as.character(speciesdata_bf$cluster)]

group_colors_bf <-  c("Sequence similarity" = "#5050FFFF", "Probabilistic taxonomic assignment"= "#CE3D32FF")  

speciesdata_bf$color <- group_colors_bf[speciesdata_bf$ID]

node <- 1:Ntip(rooted.bstree.bf)

speciesdata_bf$tip_label_new <- paste0(speciesdata_bf$Species," (", speciesdata_bf$cluster_size, ")")

bftree2024 <- (ggtree(rooted.bstree.bf, layout = 'rectangular', branch.length = TRUE) +
                 
    geom_tiplab(size =2, colour = "transparent") +
      
    geom_tippoint(aes(color = speciesdata_bf$ID[node]), size = 3) +
      
    scale_color_manual(values = group_colors_bf, name = "Identification Method") +
      
    geom_text(aes(label = speciesdata_bf$tip_label_new[node]), size =4, hjust = -0.05, fontface = "italic") +
      
    theme_tree2(legend.position = "right") +
    theme(legend.title = element_text(size = 14),legend.text = element_text(size = 12), legend.position = "top"))

bftree2024

ggsave("plots/bfphylo.png", plot = bftree2024, width = 10, height = 5, dpi = 300)


seqs_con_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 13]]
seqs_met_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 11]]
seqs_crag_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 18]]
seqs_sp_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 2]]
seqs_sub_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 17]]
seqs_tubersum_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 3]]
seqs_arct_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 14]]
seqs_vere1_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 16]]

KBIMP2024_updatedspecies <- KBIMP_updatedspecies_mos_2024 %>%
  mutate(Species = ifelse(Sample %in%names(seqs_tubersum_bf) , "Simulium tuberosum", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_vere1_bf) , "Simulium verecundum complex sp", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_con_bf) , "Simulium congareenarum", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_met_bf) , "Metacnephia borealis", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_crag_bf) , "Simulium craigi", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_sp_bf) , "Simulium sp.", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_arct_bf) , "Simulium arcticum complex sp", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_sub_bf) , "Simulium subpusillum", Species)) %>%
  filter(!Species == "unknown")


write_tsv(KBIMP2024_updatedspecies, "processed-data/KBIMP2024_updatedspecies.tsv")

rm(alighned_kbimpbf_DNA, alighned_kbimpmos_DNA, bs, bs.bf, hc, hc_bf, kbimp_bf_DNA_df,
   kbimp_bf_DNA, kbimp_mos_phydat, kbimp_mos_DNA, kbimp_bf_DNA_df, kbimp_mos_phydat, 
   KBIMP_updatedspecies_mos_2024, KBIMP2024, medoid_phy, medoid_phy_2024, medoids, 
   medoids_2024, modelTest, NJtree.kbimp.bf, NJtree.kbimp.mos, pml.tree.kbimp.bf, 
   pml.tree.kbimp.bfop, pml.tree.kbimp.mos, pml.tree.kbimp.mosop, seqs_arct_bf, 
   seqs_con_bf, seqs_crag_bf, seqs_ex_mos, seqs_ex2_mos, seqs_imp1_mos, seqs_imp2_mos, 
   seqs_imp3_mos, seqs_met_bf, seqs_punc1_mos, seqs_punc2_mos, seqs_punc3_mos, 
   seqs_punc4_mos, seqs_punc5_mos, seqs_punc6_mos, seqs_sp_bf, seqs_tubersum_bf, 
   seqs_sub_bf, seqs_vere1_bf, speciesdata, speciesdata_bf, tree_tips_bf, tree_tips_mos,
   sample_name_bf, sample_name_mos, node, groups, group_colors, group_colors_bf, 
   cluster_sizes, tree_with_bs, tree_with_bs.bf, tree_with_bsbf, dist.kbimp.bf, 
   dist.kbimp.mos, dist.medoid.bf, dist.medoid.mos, kbimp_bf_phydat)
