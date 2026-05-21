#Data required is KBIMP2025_filteredCOI.tsv from the data cleaning step and 
#the outgroup file the output files will be the updated species file and the two 
#phylogenies 


KBIMP2025 <- read_tsv(file = "../processed-data/KBIMP2025_filteredCOI.tsv")
Outgroup <- read_csv(file = "../raw-data2/Outgroup.csv")
#### PART 1 - Phylogenetic analysis and resolution of undefined Mosquito species ---- 

#filtering for Culicidae, and defining the root as a chiromidae sequence in the data set 
KBIMP2025 <- KBIMP2025 %>%
  mutate(Sample = str_replace(Sample, "^([A-Za-z]+)_(\\d+)", "\\1\\2")) %>%
  group_by(Sample) %>%
  mutate(Sample = paste0(Sample, "_", LETTERS[row_number()])) %>%
  ungroup() %>%
  bind_rows(Outgroup)

kbimp_mos_DNA_df <- KBIMP2025 %>%
  filter(Family %in% c("Culicidae", "Outgroup")) %>%
  select(Sample, Sequence) 
  
  

#preforming the alignment and converting to a PhyDat

kbimp_mos_DNA <- DNAStringSet(kbimp_mos_DNA_df$Sequence)

names(kbimp_mos_DNA) <-kbimp_mos_DNA_df$Sample

alighned_kbimpmos_DNA <- DNAStringSet(muscle::muscle(kbimp_mos_DNA))

kbimp_mos_phydat <- as.phyDat(alighned_kbimpmos_DNA, type = "DNA")

class(kbimp_mos_phydat)

length(kbimp_mos_phydat) 

###### determining medoid sequences ######

#modelTest <- modelTest(kbimp_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 
#JC is the best model because it has the lowest AIC, second lowest BIC and high loglik

dist.kbimp.mos <- dist.ml(kbimp_mos_phydat, ratio = TRUE, model = "JC69") 

#cluster sequences

hc <- hclust(dist.kbimp.mos)

groups <- cutree(hc, k = 15) 
#done multiple times to determine the optimal number of groups based on my data 

table(groups) #number of sequences in each group

cluster_sizes <- table(groups)

medoids_2025_mos <- sapply(unique(groups), function(g) {
  members <- names(groups[groups == g])
  find_medoid(dist.kbimp.mos, members)
})

medoid_phy_2025 <-kbimp_mos_phydat[medoids_2025_mos]


###### Building the tree ######

dist.medoid.mos <- dist.ml(medoid_phy_2025, ratio = TRUE, model = "JC69") 

#neighbor joining method

NJtree.kbimp.mos <- NJ(dist.medoid.mos)

plot(NJtree.kbimp.mos)

length(NJtree.kbimp.mos$tip.label) #1272 tips as expected 

# Fit the initial tree using a simple pml 

pml.tree.kbimp.mos <- pml(NJtree.kbimp.mos,kbimp_mos_phydat, k = 4, model = "GTR+I", method = "unrooted")

plot(pml.tree.kbimp.mos$tree) 

pml.tree.kbimp.mosop <- pml_bb(pml.tree.kbimp.mos, model = "GTR+I")

#bootstrapping to determine reliability

bs <- bootstrap.pml(pml.tree.kbimp.mosop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bs <- plotBS(pml.tree.kbimp.mosop$tree, bs)

#rooting the tree

rooted.bstree.mos <- root(tree_with_bs, outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.mos)

tree_with_bs <- plotBS(rooted.bstree.mos, bs)

###### Species assignment of unknown medoid groups ######

#After BLASTing these sequences I change the name of these 14 sequences and add the species name 

speciesdata <- KBIMP2025 %>%
  filter(Sample %in% c("Outgroup",
    "KGLTK0059_M_H_A" ,"KGLTK0101_M_H_A", 
    "KGLTK0275_M_H_A" ,"CBAY0169_M_H_C" ,
     "CBAY0091_M_H_B", 
    "CBAY0206_M_H_C",  "CBAY0145_M_H_A", "KGLTK0056_M_H_D",
    "KGLTK0059_M_H_B")) %>%
  select(Sample, Species) %>%
  
  mutate(ID = NA_character_) %>%

  mutate(Species = if_else(Sample == "CBAY0206_M_H_C", "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = if_else(Sample == "CBAY0145_M_H_A", "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = if_else(Sample == "KGLTK0101_M_H_A", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KGLTK0056_M_H_D", "Aedes punctor/Aedes hexodontus", Species)) %>%
  
  mutate(ID = if_else(Sample %in% c("CBAY0206_M_H_C", "CBAY0145_M_H_A", "KGLTK0101_M_H_A", "KGLTK_0056_M_H_D"), 
                      "Sequence similarity", ID)) %>%
  
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

ggsave("../plots/mosphylomoregroups.png", plot = mosquitotree2024, width = 10, height = 4, dpi = 300)


seqs_punctor_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 6]]
seqs_imp1_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 5]]
seqs_imp2_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 1]]
seqs_ex_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 7]]
seqs_un_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 9]]
seqs_punctor2_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 2]]
seqs_punctor3_mos <- alighned_kbimpmos_DNA[names(groups)[groups == 3]]

KBIMP_updatedspecies_mos <- KBIMP2025 %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punctor_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punctor2_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punctor3_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_imp1_mos) , "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_imp2_mos) , "Aedes nigripes/impiger", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_ex_mos) , "Aedes excrucians", Species)) 

#looking at the mosquito species 
  
KBIMP2025_mosquitoes <- KBIMP_updatedspecies_mos  %>%
    filter(Family == "Culicidae") %>%
    select(Species) %>%
    dplyr::count(Species, name = "countspecies")
  
length(KBIMP2025_mosquitoes$Species) #4 species of black fly are seen here but most of them are unknowns 
  
ggplot(KBIMP2025_mosquitoes, aes(y=countspecies, x=Species)) +
    geom_col(fill = "skyblue3") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

#### PART 2 - Phylogenetic analysis and resolution of undefined Blackfly species ---- 

kbimp_bf_DNA_df <- KBIMP2025 %>%
  filter(Family %in% c("Simuliidae", "Outgroup")) %>%
  select(Sample, Sequence) 

#preforming the alignment and converting to a PhyDat

kbimp_bf_DNA <- DNAStringSet(kbimp_bf_DNA_df$Sequence)

names(kbimp_bf_DNA) <-kbimp_bf_DNA_df$Sample

alighned_kbimpbf_DNA <- DNAStringSet(muscle::muscle(kbimp_bf_DNA))

kbimp_bf_phydat <- as.phyDat(alighned_kbimpbf_DNA, type = "DNA")

class(kbimp_bf_phydat)

length(kbimp_bf_phydat) 

###### determining medoid sequences ######

#modelTest <- modelTest(kbimp_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 
#JC is the best model because it has the lowest AIC, second lowest BIC and high loglik

dist.kbimp.bf <- dist.ml(kbimp_bf_phydat, ratio = TRUE, model = "JC69") 

#cluster sequences

hc <- hclust(dist.kbimp.bf)

groups <- cutree(hc, k = 30) 
#done multiple times to determine the optimal number of groups based on my data 

table(groups) #number of sequences in each group

cluster_sizes <- table(groups)

medoids_2025 <- sapply(unique(groups), function(g) {
  members <- names(groups[groups == g])
  find_medoid(dist.kbimp.bf, members)
})

medoid_phy_2025 <-kbimp_bf_phydat[medoids_2025]


###### Building the tree ######

dist.medoid.bf <- dist.ml(medoid_phy_2025, ratio = TRUE, model = "JC69") 

#neighbor joining method

NJtree.kbimp.bf <- NJ(dist.medoid.bf)

plot(NJtree.kbimp.bf)

length(NJtree.kbimp.bf$tip.label) 

# Fit the initial tree using a simple pml 

pml.tree.kbimp.bf <- pml(NJtree.kbimp.bf,kbimp_bf_phydat, k = 4, model = "GTR+I", method = "unrooted")

plot(pml.tree.kbimp.bf$tree) 

pml.tree.kbimp.bfop <- pml_bb(pml.tree.kbimp.bf, model = "GTR+I")

#bootstrapping to determine reliability

bs.bf <- bootstrap.pml(pml.tree.kbimp.bfop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bs.bf <- plotBS(pml.tree.kbimp.bfop$tree, bs.bf)

#rooting the tree

rooted.bstree.bf <- root(tree_with_bs.bf, outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.bf)

tree_with_bs <- plotBS(rooted.bstree.bf, bs.bf)

###### Species assignment of unknown medoid groups ######

#After BLASTing these sequences I change the name of these 14 sequences and add the species name 

speciesdata_bf <- KBIMP2025 %>%
  filter(Sample %in% c("CBAY0092_BF_H_A", 
                       "CBAY0095_BF_H_A" ,
                       "CBAY0147_BF_H_A" 
                       ,"CBAY0099_BF_H_E" 
                       ,"CBAY0099_BF_H_U" 
                       , "CBAY0099_BF_H_Y" 
                       , "CBAY0103_BF_H_B" 
                       , "CBAY0104_BF_H_A" 
                       , "CBAY0278_BF_H_C" 
                       , "CBAY0104_BF_H_D" 
                       ,"CBAY0157_BF_H_B" 
                       , "CBAY0132_BF_H_A" 
                       , "CBAY0166_BF_H_C" 
                       , "CBAY0201_BF_H_B" 
                       ,"CBAY0247_BF_H_A" 
                       ,"CBAY0267_BF_H_D" 
                       , "CBAY0277_BF_H_B" 
                       , "CBAY0306_BF_H_B" 
                       ,"KGLTK0113_BF_H_E"
                       , "KGLTK0151_BF_H_A"
                       , "KGLTK0111_BF_H_A"
                       , "KGLTK0103_BF_H_A"
                       ,"KGLTK0107_BF_H_A"
                       ,"KGLTK0107_BF_H_C"
                       , "KGLTK0112_BF_H_A"
                       ,"KGLTK0137_BF_H_A"
                       ,"KGLTK0137_BF_H_B"
                       , "KGLTK0137_BF_H_E"
                       , "KGLTK0145_BF_H_A"
                       , "Outgroup"        
  )) %>%
  select(Sample, Species) %>%
  
  mutate(ID = NA_character_) %>%
  
  mutate(Species = if_else(Sample == "CBAY0095_BF_H_A", "Simulidae sp.", Species)) %>%
  mutate(Species = if_else(Sample == "KGLTK0107_BF_H_A", "Simulium verecundum complex sp", Species)) %>%
  mutate(Species = if_else(Sample == "KGLTK0107_BF_H_C", "Simulium verecundum complex sp", Species)) %>%
  mutate(Species = if_else(Sample == "KGLTK0111_BF_H_A", "Metacnephia bilineata", Species)) %>%
  mutate(Species = if_else(Sample == "CBAY0132_BF_H_A", "Simulium tuberosum", Species)) %>%
  mutate(Species = if_else(Sample == "CBAY0099_BF_H_U", "Metacnephia borealis", Species)) %>%
  mutate(Species = if_else(Sample == "CBAY0099_BF_H_Y", "Metacnephia borealis", Species)) %>%
  mutate(Species = if_else(Sample == "CBAY0145_BF_H_A", "Simulium arcticum complex sp", Species)) %>%
  
  mutate(ID = if_else(Sample %in% c("CBAY0095_BF_H_A", "KGLTK0107_BF_H_A",
                                    "KGLTK0107_BF_H_C","KGLTK0111_BF_H_A",
                                    "CBAY0132_BF_H_A", "CBAY0099_BF_H_U",
                                    "CBAY0145_BF_H_A", "CBAY0099_BF_H_Y"), 
                      "Sequence similarity", ID)) %>%
  
  mutate(ID = if_else(is.na(ID), "Probabilistic taxonomic assignment", ID))


##### Finalized Figure ######

# match trait data and tree data

tree_tips_bf <- rooted.bstree.bf$tip.label 

sample_name_bf <- speciesdata_bf$Sample 

# checking if all species names are present in the tree

all(tree_tips_bf %in% sample_name_bf) 

all(sample_name_bf %in% tree_tips_bf)

setdiff(tree_tips_bf, sample_name_bf)

setdiff(sample_name_bf, tree_tips_bf)

speciesdata_bf <- speciesdata_bf[match(tree_tips_bf, speciesdata_bf$Sample), ] %>%
  rownames_to_column(var = "cluster")

speciesdata_bf$cluster_size <- cluster_sizes[as.character(speciesdata_bf$cluster)]

group_colors <-  c("Sequence similarity" = "#5050FFFF", "Probabilistic taxonomic assignment" = "#CE3D32FF")

speciesdata_bf$color <- group_colors[speciesdata_bf$ID]

#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.bf)

speciesdata_bf$tip_label_new <- paste0(speciesdata_bf$Species," (", speciesdata_bf$cluster_size, ")")

#produce the tree

bftree2024 <- (
  ggtree(rooted.bstree.bf, layout = "rectangular", branch.length = TRUE) +
    
    geom_tippoint(aes(color = speciesdata_bf$ID[node]), size = 3) +
    
    scale_color_manual(values = group_colors, name = "Identification Method") +
    
    geom_text(aes(label = speciesdata_bf$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
    
    theme(legend.title = element_text(size = 14),      
          legend.text = element_text(size = 12), legend.position = "top")
)

plot(bftree2024)

#save the tree

ggsave("../plots/mosphylomoregroups.png", plot = mosquitotree2024, width = 10, height = 4, dpi = 300)


seqs_tubersum_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 12]]
seqs_vere1_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 23]]
seqs_vere2_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 24]]
seqs_sub_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 2]]
seqs_meta1_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 5]]
seqs_meta2_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 6]]
seqs_meta3_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 3]]
seqs_bil_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 21]]
seqs_arct_bf <- alighned_kbimpbf_DNA[names(groups)[groups == 29]]



KBIMP2025_updatedspecies <- KBIMP_updatedspecies_mos %>%
  mutate(Species = ifelse(Sample %in%names(seqs_tubersum_bf) , "Simulium tuberosum", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_vere1_bf) , "Simulium verecundum complex sp", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_vere2_bf) , "Simulium verecundum complex sp", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_meta1_bf) , "Metacnephia borealis", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_meta2_bf) , "Metacnephia borealis", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_meta3_bf) , "Metacnephia borealis", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_sub_bf) , "Simulium sp.", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_arct_bf) , "Simulium arcticum complex sp", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_bil_bf) , "Metacnephia bilineata", Species)) %>%
  mutate(Species = str_replace(Species, "Simulium tuberosum complex", "Simulium tuberosum")) %>%
  filter(!Species == "unknown")

#looking at the mosquito species 

KBIMP2025_mosquitoes <- KBIMP_updatedspecies  %>%
  filter(Family == "Simuliidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP2025_mosquitoes$Species) #4 species of black fly are seen here but most of them are unknowns 

ggplot(KBIMP2025_mosquitoes, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

write_tsv(KBIMP2025_updatedspecies, "../processed-data/KBIMP2025_updatedspecies.tsv")
