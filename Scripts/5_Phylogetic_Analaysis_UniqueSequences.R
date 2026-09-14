
#Data required for analysis is filtered data from 2024 and 2025 (KBIMP2024_filteredCOI and KBIMP2025_filteredCOI) 
#KBIMP 2024 sample data to fix up the 2024 species data 




#### PART 1 - Loading packages, Openeing files, Running functions ----

##### Packages #####

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
library(gt)

##### Files #####

KBIMP2024 <- read_tsv(file = "processed-data/KBIMP2024_filteredCOI.tsv")
KBIMP2025 <- read_tsv(file = "processed-data/KBIMP2025_filteredCOI.tsv")
kbimp2024_sampledata_clean <- read_csv( 
          files = "processed-data/kbimp2024_sampledata_clean.csv")
Outgroup <- read_csv(file = "raw-data2/Outgroup.csv")


##### Packages #####

# function to find medoid sequence in a cluster
find_medoid <- function(distmat, members) {
  mat <- as.matrix(distmat)
  members <- intersect(members, rownames(mat)) # only keep matching names
  submat <- mat[members, members, drop = FALSE]
  sums <- rowSums(submat)
  medoid <- names(which.min(sums))
  return(medoid)
}

#binary trait matrix

get_binary_trait_matrix <- function(seq_data) {
  
  mosmetadata <- kbimp_mosmetadata %>%
    filter(SampleID %in% names(seq_data)) %>%
    select(SampleID, ExactSite) %>%
    column_to_rownames(var = "SampleID")
  
  binary_matrix_trait <- model.matrix(~ ExactSite - 1, data = mosmetadata)
  
  binary_matrix_trait <- binary_matrix_trait %>% as.data.frame() %>%
    rownames_to_column(var = "SampleID") 
  
  return(binary_matrix_trait)
  
}


#### PART 2 - Combining data from 2024 and 2025 ----

#preparing data from 2024 to be in the same format 

KBIMP2024 <- KBIMP2024 %>%
  left_join(kbimp2024_sampledata_clean, 
            join_by(Sample == SampleID)) %>%
  mutate(Sample = paste0(FieldID, case_when(
                             Family == "Simuliidae" ~ "_BF",
                             Family == "Culicidae" ~ "_M", TRUE ~ ""))) %>%
  select(Sample, Read_Count, Kingdom, Phylum, Class, Order, Family, Genus, Species, Probability_Genus, Probability_Species, Sequence) %>%
  group_by(Sample) %>%
  mutate(Sample = paste0(Sample, "_", LETTERS[row_number()])) %>%
  ungroup()

#writing this file to be used in the sample count 

write_tsv(KBIMP2024, "processed-data/KBIMP2024_clean.tsv")

#preparing 2025 data for combination 

KBIMP2025 <- KBIMP2025 %>%
  mutate(
    Sample = Sample %>%
      str_remove("_KFS$") %>%   
      str_remove("_H$") %>%   
      str_replace("^([A-Za-z]+)_(\\d+)", "\\1\\2")) %>%

  group_by(Sample) %>%
  mutate(Sample = paste0(Sample, "_", LETTERS[row_number()])) %>%
  ungroup() %>%
    select(Sample, Read_Count, Kingdom, Phylum, Class, Order, Family, Genus, Species, Probability_Genus, Probability_Species, Sequence) 

#combining 

KBIMP <- bind_rows(KBIMP2024, KBIMP2025, Outgroup)

#### tree both years mosquitoes ---- 

kbimp_mos_DNA_df <- KBIMP %>%
  filter(Family %in% c("Culicidae", "Outgroup")) %>%
  select(Sample, Sequence) 

##### aligning and preparing phydat #####

kbimp_mos_DNA <- DNAStringSet(kbimp_mos_DNA_df$Sequence)

names(kbimp_mos_DNA) <- kbimp_mos_DNA_df$Sample

alighned_kbimpmos_DNA <- DNAStringSet(muscle::muscle(kbimp_mos_DNA))

#browseSeqs(aligned_kbimpmos_DNA)

unique_seqsmos <- unique(alighned_kbimpmos_DNA)

writeXStringSet(unique_seqsmos,
                filepath = "processed-data/uniquemosseqforbold.fasta",
                format = "fasta")

kbimp_mos_phydat <- as.phyDat(unique_seqsmos, type = "DNA")
class(kbimp_mos_phydat) # is a "phyDat" object
length(kbimp_mos_phydat) # 113 has the  species as seen before

modelTest <- modelTest(kbimp_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) 

###### Building the tree ######

dist.mos <- dist.ml(kbimp_mos_phydat, ratio = TRUE) 

#neighbor joining method

NJtree.kbimp.mos <- NJ(dist.mos)

plot(NJtree.kbimp.mos)

length(NJtree.kbimp.mos$tip.label) #68 tips as expected 

# Fit the initial tree using a simple pml 

pml.tree.kbimp.mos <- pml(NJtree.kbimp.mos, kbimp_mos_phydat, k = 4, model = "GTR", method = "unrooted")

plot(pml.tree.kbimp.mos$tree) 

pml.tree.kbimp.mosop <- pml_bb(pml.tree.kbimp.mos, 
                               model = "GTR")

#bootstrapping to determine reliability

bs <- bootstrap.pml(pml.tree.kbimp.mosop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bs <- plotBS(pml.tree.kbimp.mos$tree, bs)

#rooting the tree

rooted.bstree.mos <- root(tree_with_bs, outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.mos)

tree_with_bs <- plotBS(rooted.bstree.mos, bs)

##### preparing labels for finalzed tree #####

# brining in BOLD data for tree 

BOLDID_mos <- read.csv(file = "processed-data/BOLDID_aedes.csv")

BOLDID_mos2 <- BOLDID_mos  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[(BOLD:.*?)\\]",
          remove = FALSE) %>% 
  
  filter(!is.na(BOLDID)) %>%
  
  filter(ID. >= 97) %>% 
  
  group_by(Query.ID) %>%
  
  slice_max(ID., n = 1, with_ties = FALSE) %>%
  
  ungroup() 

# trait data and the tree have the same order of species 

tree_tips_mos <- rooted.bstree.mos$tip.label

sample_name_mos <- BOLDID_mos2$Query.ID 

# checking if all species names are present in the tree

all(tree_tips_mos %in% sample_name_mos)

all(sample_name_mos %in% tree_tips_mos) 

setdiff(tree_tips_mos, sample_name_mos)

setdiff(sample_name_mos, tree_tips_mos)

#getting trait data set up as tip labels 

BOLDID_mos2 <- BOLDID_mos2[match(tree_tips_mos, BOLDID_mos2$Query.ID), ] 

BOLDID_mos2$tip_label_new <- paste0(BOLDID_mos2$BOLDID," (", BOLDID_mos2$ID., ")")

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.mos$node.label)

internal_nodes <- (Ntip(rooted.bstree.mos)+1):(Ntip(rooted.bstree.mos)+Nnode(rooted.bstree.mos))

bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  
  mutate(bootstrap = (bootstrap *100)) %>%
  
  filter(bootstrap >= 60) 

#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.mos)

##### Produce finalized tree #####

mosquitotree2024 <- (ggtree(rooted.bstree.mos, layout = "rectangular", branch.length = TRUE) +
                       
                geom_text(aes(label = ifelse(label =="Outgroup",
                                                  "Outgroup", 
                          BOLDID_mos2$tip_label_new[node])),
                          hjust = -0.05, size = 4) +
        
        geom_strip('KGLTK0008_M_J', 'KGLTK0018_M_C', barsize=2, 
                   color='skyblue3', hjust = 1.2,  fontface = "italic",
                   label="Culiseta alaskaensis", offset.text=.1, fontsize =5) +
        geom_strip('KGLTK0101_M_A', 'KGLTK0018_M_A', barsize=2, 
                   color='pink2', hjust = 1.3, fontface = "italic",
                   label="Aedes excrucians", offset.text=.1, fontsize =5) +
        geom_strip('CBAY0036_M_I', 'CBAY0206_M_C', barsize=2, 
                   color='gold', hjust = 1.2, fontface = "italic",
                   label="Aedes nigripes/impiger", 
                   offset.text=.1, fontsize =5) +
        geom_strip('CBAY0300_M_A', 'CBAY0108_M_A', barsize=2, 
                   color='green3', hjust = 1.3, fontface = "italic",
                   label="Aedes communis", offset.text=.1, fontsize =5) +
        geom_strip('KGLTK0008_M_E', 'CBAY0213_M_C', barsize=2, 
                   color='darkblue', hjust = 1, fontface = "italic",
                   label="Aedes punctor/Aedes hexodontus", 
                   offset.text=.1, fontsize =5) +
        geom_strip('KGLTK0008_M_L', 'KGLTK0008_M_L', barsize=2, 
                   color='blueviolet', hjust = 1.3, fontface = "italic",
                   label="Culiseta inornata", offset.text=.1, fontsize =5) +
                       
      theme(legend.title = element_text(size = 14), 
          legend.text = element_text(size = 12), 
          legend.position = "top")) %<+% 
  
  bs_tibble + geom_label2(aes(label = bootstrap), 
                          hjust = 0.7, size = 3, 
                          color = "red", fill = "white") 

ggsave("plots/mosallseqtree.png", plot = mosquitotree2024, width = 14, height = 16, dpi = 300)


#### both years black flies just Simulium genus ----

kbimp_sim_DNA_df <- KBIMP %>%
  filter(Genus %in% c("Simulium", "Outgroup")) %>%
  select(Sample, Sequence)  

##### aligning and preparing phydat #####

kbimp_sim_DNA <- DNAStringSet(kbimp_sim_DNA_df$Sequence)

names(kbimp_sim_DNA) <- kbimp_sim_DNA_df$Sample

alighned_kbimpsim_DNA <- DNAStringSet(muscle::muscle(kbimp_sim_DNA))

#BrowseSeqs(alighned_kbimpbf_DNA)

unique_seqssim <- unique(alighned_kbimpsim_DNA)

writeXStringSet(unique_seqssim,
                filepath = "processed-data/unique_simseq_forbold.fasta",
                format = "fasta")

kbimp_sim_phydat <- as.phyDat(unique_seqssim, type = "DNA")
class(kbimp_sim_phydat) # is a "phyDat" object
length(kbimp_sim_phydat) # 83 unique seq

#maybe should do this once at the begginig on the other file 
modelTest <- modelTest(kbimp_sim_phydat, model=c("JC", "JC69", "F81", "K80", "HKY", "SYM", "GTR", "TN93", "GG95")) 

###### Building the tree ######

#create a new dist matrix

dist.medoid.sim <- dist.ml(kbimp_sim_phydat, ratio = TRUE, model = "JC69") 

# creating a tree using the neighbor joining method

NJtree.kbimp.sim <- NJ(dist.medoid.sim)

plot(NJtree.kbimp.sim)

length(NJtree.kbimp.sim$tip.label) 

# Fit the initial tree using a simple pml

pml.tree.kbimp.sim <- pml(NJtree.kbimp.sim, kbimp_sim_phydat, k = 4, model = "GTR+I", method = "unrooted")

plot(pml.tree.kbimp.sim$tree) 

pml.tree.kbimp.simop <- pml_bb(pml.tree.kbimp.sim, model = "GTR+I")

#bootstrapping analysis for tree

bs.sim <- bootstrap.pml(pml.tree.kbimp.simop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bssim <- plotBS(pml.tree.kbimp.sim$tree, bs.sim)

#rooting the 

rooted.bstree.sim <- root(tree_with_bssim, 
                          outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.sim)

tree_with_bs.sim <- plotBS(rooted.bstree.sim, bs.sim)

##### preparing labels for finalzed tree #####

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.sim$node.label)

internal_nodes <- (Ntip(rooted.bstree.sim)+1):(Ntip(rooted.bstree.sim)+Nnode(rooted.bstree.sim))

bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  mutate(bootstrap = (bootstrap *100)) %>%
  filter(bootstrap >= 60) 

#bringing in bold data for final tree 

BOLDID_sim <- read.csv(file = "processed-data/BOLDID_simuulidae.csv")

BOLDID_sim2 <- BOLDID_sim  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[(BOLD:.*?)\\]",
          remove = FALSE) %>% 
  
  filter(!is.na(BOLDID)) %>%
  
  filter(ID. >= 94) %>%
  
  group_by(Query.ID) %>%
  
  slice_max(ID., n = 1, with_ties = FALSE) %>%
  
  ungroup() 

# trait data and the tree have the same order of species 

tree_tips_sim <- rooted.bstree.sim$tip.label

sample_name_sim <- BOLDID_sim2$Query.ID 

# checking if all species names are present in the tree

all(tree_tips_sim %in% sample_name_sim)

all(sample_name_sim %in% tree_tips_sim) 

setdiff(tree_tips_sim, sample_name_sim)

setdiff(sample_name_sim, tree_tips_sim)

#getting trait data set up as tip labels 

BOLDID_sim2 <- BOLDID_sim2[match(tree_tips_sim, BOLDID_sim2$Query.ID), ] 

BOLDID_sim2$tip_label_new <- paste0(BOLDID_sim2$BOLDID," (", BOLDID_sim2$ID., ")")


#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.sim)

##### Produce finalized tree #####

simtreespecies2024 <- (ggtree(rooted.bstree.sim, layout = "rectangular", branch.length = TRUE) +
                        
              geom_text(aes(label = BOLDID_sim2$tip_label_new[node]),
               hjust = -0.05, size =4, fontface = "italic") +
                        
             geom_strip('KGLTK0049_BF_C', 'KGLTK0049_BF_B', 
              barsize=2, color='skyblue3',  
              fontface = "italic", offset = -0.01,
              label= "Simulium arcticum complex", 
              offset.text=0.003, fontsize =5) + 
            geom_strip('KGLTK0036_BF_C', 'KGLTK0146_BF_A', 
              barsize=2, color='pink2',  fontface = "italic",
              label= "Simulium malyschevi",  offset = -0.01,
              offset.text=.003, fontsize =5) + 
              geom_strip('KGLTK0012_BF_B', 'CBAY0308_BF_C', 
                        barsize=2, color='gold',  
                         fontface = "italic",offset = -0.01,
                         label= "Simulium decimatum", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('KGLTK0137_BF_C', 'KGLTK0137_BF_B', 
                        barsize=2, color='darkblue',  
                         fontface = "italic", offset = -0.01,
                         label= "Simulium murmanum", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('KGLTK0019_BF_C', 'KGLTK0034_BF_B', 
                         barsize=2, color='green3',  offset = -0.01,
                         fontface = "italic",
                         label= "Simulium venustum complex", 
                         offset.text=.002, fontsize =5) + 
              geom_strip('CBAY0078_BF_A', 'CBAY0267_BF_D', 
                         barsize= 2, color='darkolivegreen1',  
                         fontface = "italic", offset = 0.004,
                         label= "Simulium noelleri", 
                         offset.text=.002, fontsize =5) + 
              geom_strip('CBAY0062_BF_B', 'CBAY0062_BF_B', 
                         barsize=2, color='orange',  
                         fontface = "italic",
                         label= "Simulium decorum", 
                         offset.text=-0.01, fontsize =5) + 
              geom_strip('KGLTK0107_BF_A', 'KGLTK0107_BF_A', 
                         barsize=2, color='red3',  
                         fontface = "italic", offset = -0.01,
                         label= "Simulium verecundum complex", 
                         offset.text=.0006, fontsize =5) + 
              geom_strip('CBAY0169_BF_A', 'CBAY0106_BF_E', 
                         barsize=2, color='forestgreen',  
                         fontface = "italic", offset = -0.01,
                         label= "Simulium tuberosum", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('CBAY0275_BF_C', 'CBAY0275_BF_G', 
                         barsize=2, color='maroon',  
                         fontface = "italic", offset = 0.01,
                         label= "Simulium baffinense", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('KGLTK0137_BF_E', 'KGLTK0137_BF_D', 
                         barsize=2, color='brown2',   
                         fontface = "italic", offset = -0.001,
                         label= "Simulium silvestre", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('CBAY0056_BF_B', 'CBAY0277_BF_C', 
                         barsize=2, color='seagreen2',   
                         fontface = "italic", offset = -0.01,
                         label= "Simulium subpusillum", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('CBAY0157_BF_C', 'CBAY0157_BF_B', 
                         barsize=2, color='cyan',  
                         fontface = "italic", offset = -0.01,
                         label= "Simulium craigi", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('KGLTK0112_BF_B', 'KGLTK0015_BF_B', 
                         barsize=2, color='blueviolet',  
                         fontface = "italic", offset = -0.01,
                         label= "Simulium congareenarum", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('KGLTK0103_BF_A', 'KGLTK0103_BF_A', 
                         barsize=2, color='darkgoldenrod3',  
                         fontface = "italic",
                         label= "Simulium excisum", 
                         offset.text=-.001, fontsize =5) + 
              geom_strip('CBAY0222_BF_A', 'CBAY0306_BF_B', 
                         barsize=2, color='deeppink2',   
                         fontface = "italic", offset = -0.01,
                         label= "Simulium vittatum", 
                         offset.text=.003, fontsize =5) + 
              geom_strip('Outgroup', 'Outgroup', 
                         barsize=2, color='black',   
                         fontface = "italic",
                         label= "Outgroup", 
                         offset.text=.05, fontsize =5) + 
                        
      theme(legend.title = element_text(size = 14), 
      legend.text = element_text(size = 12), 
      legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, 
              color = "red", fill = "white") 

simtreespecies2024 

ggsave("plots/Simulium_treespecies.png", plot = simtreespecies2024, width = 15, height = 17, dpi = 300)


#### both years black flies not Simulium genus ----

kbimp_meta_DNA_df <- KBIMP %>%
  filter(Genus %in% c("Metacnephia", "Outgroup")) %>%
  select(Sample, Sequence)  

##### aligning and preparing phydat #####

kbimp_meta_DNA <- DNAStringSet(kbimp_meta_DNA_df$Sequence)

names(kbimp_meta_DNA) <- kbimp_meta_DNA_df$Sample

alighned_kbimpmeta_DNA <- DNAStringSet(muscle::muscle(kbimp_meta_DNA))

#BrowseSeqs(alighned_kbimpbf_DNA)

unique_seqsmeta <- unique(alighned_kbimpmeta_DNA)

writeXStringSet(unique_seqsmeta,
                filepath = "processed-data/unique_metaseq_forbold.fasta",
                format = "fasta")

kbimp_meta_phydat <- as.phyDat(unique_seqsmeta, type = "DNA")
class(kbimp_meta_phydat) # is a "phyDat" object
length(kbimp_meta_phydat) # 83 unique seq

###### Building the tree ######


#bringing in bold data for final tree 

BOLDID_meta <- read.csv(file = "processed-data/BOLDID_meta.csv")

BOLDID_meta2 <- BOLDID_meta  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[(BOLD:.*?)\\]",
          remove = FALSE) %>% 
  
  filter(!is.na(BOLDID)) %>%
  
  filter(ID. >= 94) %>%
  
  group_by(Query.ID) %>%
  
  slice_max(ID., n = 1, with_ties = FALSE) %>%
  
  ungroup() 


#create a new dist matrix

dist.medoid.meta <- dist.ml(kbimp_meta_phydat, ratio = TRUE, model = "JC69") 

# creating a tree using the neighbor joining method

NJtree.kbimp.meta <- NJ(dist.medoid.meta)

plot(NJtree.kbimp.meta)

length(NJtree.kbimp.meta$tip.label) 

# Fit the initial tree using a simple pml

pml.tree.kbimp.meta <- pml(NJtree.kbimp.meta, kbimp_meta_phydat, k = 4, model = "GTR+I", method = "unrooted")

plot(pml.tree.kbimp.meta$tree) 

pml.tree.kbimp.metaop <- pml_bb(pml.tree.kbimp.meta, model = "GTR+I")

#bootstrapping analysis for tree

bs.meta <- bootstrap.pml(pml.tree.kbimp.metaop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bsmeta <- plotBS(pml.tree.kbimp.meta$tree, bs.meta)

#rooting the 

rooted.bstree.meta <- root(tree_with_bsmeta, outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.meta)

tree_with_bs.meta <- plotBS(rooted.bstree.meta, bs.meta)

##### preparing labels for finalzed tree #####

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.meta$node.label)

internal_nodes <- (Ntip(rooted.bstree.meta)+1):(Ntip(rooted.bstree.meta)+Nnode(rooted.bstree.meta))

bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  mutate(bootstrap = (bootstrap *100)) %>%
  filter(bootstrap >= 60) 


#trait data and the tree have the same order of species 

tree_tips_meta <- rooted.bstree.meta$tip.label

sample_name_meta <- BOLDID_meta2$Query.ID 

# checking if all species names are present in the tree

all(tree_tips_meta %in% sample_name_meta)

all(sample_name_meta %in% tree_tips_meta) 

setdiff(tree_tips_meta, sample_name_meta)

setdiff(sample_name_meta, tree_tips_meta)

#getting trait data set up as tip labels 

BOLDID_meta2 <- BOLDID_meta2[match(tree_tips_meta, BOLDID_meta2$Query.ID), ] 

BOLDID_meta2$tip_label_new <- paste0(BOLDID_meta2$BOLDID," (", BOLDID_meta2$ID., ")")


#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.meta)

##### Produce finalized tree #####

metatreespecies2024 <- (ggtree(rooted.bstree.meta, layout = "rectangular",
                               branch.length = TRUE) +
                          
                          geom_text(aes(label = BOLDID_meta2$tip_label_new[node]), 
                                    hjust = -0.05, size =4, fontface = "italic") +
                          
                          geom_strip('CBAY0073_BF_J', 'CBAY0071_BF_A', 
                                     barsize=2, color='skyblue2',  
                                     fontface = "italic", offset = -0.05,
                                     label= "Metacnephia borealis", 
                                     offset.text=.003, fontsize =5) +
                          geom_strip('CBAY0264_BF_D', 'CBAY0294_BF_C', 
                                     barsize=2, color='skyblue2',  
                                     fontface = "italic", offset = -0.01,
                                     label= "Metacnephia borealis", 
                                     offset.text=.003, fontsize =5) +
                          geom_strip('KGLTK0111_BF_A', 'KGLTK0113_BF_D', 
                                     barsize=2, color='gold',  
                                     fontface = "italic", offset = -0.007,
                                     label= "Metacnephia bilineata", 
                                     offset.text=.003, fontsize =5) +
                          geom_strip('Outgroup', 'Outgroup', 
                                     barsize=2, color='black',   
                                     fontface = "italic", 
                                     label= "Outgroup", offset = 0.01,
                                     fontsize =5) +
                          
                          theme(legend.title = element_text(size = 14), 
                                legend.text = element_text(size = 12), 
                                legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, 
              color = "red", fill = "white") 

metatreespecies2024 


ggsave("plots/meta_treespecies.png", plot = metatreespecies2024, width = 18, height = 5, dpi = 300)


##### species assignments -----

BOLDIDspecies <- read_csv(file = "processed-data/BOLDIDspecies.csv")

# Combine BOLD tables
BOLDresults <- bind_rows(BOLDID_mos2, BOLDID_sim2, BOLDID_meta2)

BOLDresults <- BOLDresults %>%
  filter(!Query.ID == "Outgroup") %>%
  filter(!is.na(Query.ID)) %>%
  select(Query.ID, BOLDID)

KBIMP_updatedspecies <- KBIMP %>%
  filter(Family %in% c("Culicidae", "Simuliidae")) %>%
  mutate(Species = str_replace(Species, "Aedes punctor", "Aedes punctor/Aedes hexodontus")) %>%
  left_join(BOLDresults, join_by("Sample" =="Query.ID")) %>%
  group_by(Sequence) %>%
  mutate(BOLDID = ifelse(any(!is.na(BOLDID)), na.omit(BOLDID)[1], NA)) %>%
  ungroup() %>%
  left_join(BOLDIDspecies, relationship = "many-to-many", 
            join_by(BOLDID == BIN)) %>%
  mutate(update_flag = case_when(
      Species == "unknown" | is.na(Species) ~ "Sequence Simularity",
      TRUE ~ "Probabolistic")) %>%
    mutate(Species = coalesce(Species_BOLDID, Species)) %>%
    select(-Species_BOLDID) %>%
    filter(!Species == "unknown") %>%
  select(Sample, Species, Genus, Family, update_flag, Sequence) 
 

write_tsv(KBIMP_updatedspecies, "processed-data/KBIMP_updatedspecies.tsv")


BOLDinvest <- bind_rows(BOLDID_mos, BOLDID_sim, BOLDID_meta)

BOLDinvest <- BOLDinvest %>%
  right_join(KBIMP_updatedspecies, join_by("Query.ID" == "Sample")) %>%
filter(Species.x == Species.y) %>%
  group_by(Query.ID) %>%
  slice_max(ID., n = 1, with_ties = FALSE) 

BOLDIDspecies_table <- BOLDIDspecies %>%
  mutate(`Percentage of species assignments` = 
           ifelse(total_num_species == 0, NA, 
                  num_species_repo / total_num_species *100)) %>%
  mutate(`Species assignment:Total sequences in BIN` = paste0(num_species_repo, ":", total_observations)) %>%
  dplyr::rename(Species = Species_BOLDID) %>%
  select(BIN, Species, 
         `Percentage of species assignments`, 
         `Species assignment:Total sequences in BIN`) %>% 
  arrange(Species) %>%
  gt() %>%
  fmt_missing(columns = everything(),
              missing_text = "") %>%
  fmt_number(columns = 2:last_col(), 
             decimals = 2) %>%
  cols_align(align = "center",
    columns = c(`Percentage of species assignments`, 
                `Species assignment:Total sequences in BIN`)) %>%
  tab_stubhead(label = md("")) %>%
  tab_header(title = md("")) %>%
  tab_footnote(footnote = md("")) %>%
  tab_options(
    table.border.top.width = px(0),
    table.border.bottom.width = px(0),
    column_labels.border.top.width = px(0),
    column_labels.border.bottom.width = px(0),
    table_body.hlines.width = px(0),
    table_body.vlines.width = px(0),
    row_group.border.top.width = px(0),
    row_group.border.bottom.width = px(0),
    stub.border.style = "none") %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = list(cells_column_spanners(), 
                             cells_row_groups())) %>%
  tab_style(style = cell_borders(sides = c("bottom"),
                                 color = "black", weight = px(2)),
            locations = cells_title()) %>%
  tab_style(style = cell_borders(sides = c("top"),
                                 color = "black", weight = px(2)),
            locations = cells_footnotes()) %>%
  tab_style(style = cell_borders(sides = c("bottom"),
                                 color = "black", weight = px(2)),
            locations = list(cells_column_labels(), 
                             cells_stubhead())) %>%
  tab_style(style = cell_borders(
    sides = c("top", "bottom"),
    color = "black", weight = px(2)),
    locations = cells_row_groups()) %>%
  tab_options(data_row.padding = px(5)) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(Species)))
  
gtsave(data = BOLDIDspecies_table , 
       filename = "plots/BOLDIDspecies_table.png")


#### Analysis of genetic divergences ----

##### black flies #####

#filter out just Simulium

bf_furtherinvestigation <- KBIMP_updatedspecies %>%
  filter(Genus == "Simulium") %>%
  select(Species, Sequence) 

#create DNA stringset

DNA_furtherinvestigation_bf <- DNAStringSet(bf_furtherinvestigation$Sequence)

names(DNA_furtherinvestigation_bf) <- bf_furtherinvestigation$Species

alighned_DNA_furtherinvestigation_bf <- DNAStringSet(muscle::muscle(DNA_furtherinvestigation_bf))

unique_seqs_bf <- unique(alighned_DNA_furtherinvestigation_bf)

kbimp_bf_phydat <- as.phyDat(unique_seqs_bf, type = "DNA")

#calculate distance matrix and assosiate this with the species names

dist.bf <- dist.ml(kbimp_bf_phydat, model = "F81", ratio = TRUE) 

species_bf <- names(unique_seqs_bf)

dist_mat <- as.matrix(dist.bf)

#find the average distance for each species pairing 

averagedist_bf <- data.frame()
unique_species_bf <- unique(species_bf)

# Start j at 'i' instead of 'i+1' to include within-species
for (i in 1:length(unique_species_bf)) {
  for (j in i:length(unique_species_bf)) {
    
    sp1 <- unique_species_bf[i]
    sp2 <- unique_species_bf[j]
    
    idx1 <- which(species_bf == sp1)
    idx2 <- which(species_bf == sp2)
    
    # Extract the sub-matrix for these two groups
    vals <- dist_mat[idx1, idx2]
    
    # CRITICAL: If comparing the same species, remove the diagonal (self-comparisons)
    if (sp1 == sp2) {
      # If there's only 1 individual in the species, mean is NA (can't compare to others)
      if (length(idx1) == 1) {
        mean_val <- NA
      } else {
        # Convert to matrix and extract only the off-diagonal elements
        vals_mat <- as.matrix(vals)
        mean_val <- mean(vals_mat[row(vals_mat) != col(vals_mat)], na.rm = TRUE)
      }
    } else {
      # For different species, just take the standard mean
      mean_val <- mean(as.vector(vals), na.rm = TRUE)
    }
    
    averagedist_bf <- rbind(averagedist_bf, data.frame(
      Species1 = sp1,
      Species2 = sp2,
      MeanDistance = mean_val
    ))
  }
}

#create final table to export

bf_dis_matrix <- averagedist_bf %>%
  pivot_wider(names_from = Species2, values_from = MeanDistance) %>%
  column_to_rownames("Species1") %>%
  as_tibble(rownames = "Species") %>%
  gt() %>%
  fmt_missing(columns = everything(),
              missing_text = "") %>%
  fmt_number(columns = 2:last_col(), 
             decimals = 5) 

gtsave(data = bf_dis_matrix , 
       filename = "plots/bf_dis_matrix.png", vwidth = 2300, vheight = 1500)

##### mosquitoes #####

#filter out mosquitoes from data 

mos_furtherinvestigation <- KBIMP_updatedspecies %>%
  filter(Family == "Culicidae") %>%
  select(Species, Sequence) 

#create DNa string set/ phydat

DNA_furtherinvestigation <- DNAStringSet(mos_furtherinvestigation$Sequence)

names(DNA_furtherinvestigation) <- mos_furtherinvestigation$Species

alighned_DNA_furtherinvestigation <- DNAStringSet(muscle::muscle(DNA_furtherinvestigation))

unique_seqs <- unique(alighned_DNA_furtherinvestigation)

kbimp_mos_phydat <- as.phyDat(unique_seqs, type = "DNA")

#calculate distances 

dist.mos <- dist.ml(kbimp_mos_phydat, ratio = TRUE) 

#assosiate distances with species 

species <- names(unique_seqs)

dist_mat <- as.matrix(dist.mos)

#get the average values 

results <- data.frame()

unique_species <- unique(species)


# Start j at 'i' instead of 'i+1' to include within-species
for (i in 1:length(unique_species)) {
  for (j in i:length(unique_species)) {
    
    sp1 <- unique_species[i]
    sp2 <- unique_species[j]
    
    idx1 <- which(species == sp1)
    idx2 <- which(species == sp2)
    
    # Extract the sub-matrix for these two groups
    vals <- dist_mat[idx1, idx2]
    
    # CRITICAL: If comparing the same species, remove the diagonal (self-comparisons)
    if (sp1 == sp2) {
      # If there's only 1 individual in the species, mean is NA (can't compare to others)
      if (length(idx1) == 1) {
        mean_val <- NA
      } else {
        # Convert to matrix and extract only the off-diagonal elements
        vals_mat <- as.matrix(vals)
        mean_val <- mean(vals_mat[row(vals_mat) != col(vals_mat)], na.rm = TRUE)
      }
    } else {
      # For different species, just take the standard mean
      mean_val <- mean(as.vector(vals), na.rm = TRUE)
    }
    
    results <- rbind(results, data.frame(
      Species1 = sp1,
      Species2 = sp2,
      MeanDistance = mean_val
    ))
  }
}


mosquito_dis_matrix <- results %>%
  pivot_wider(names_from = Species2, values_from = MeanDistance) %>%
  column_to_rownames("Species1") %>%
  as_tibble(rownames = "Species") %>%
  gt() %>%
  fmt_missing(columns = everything(),
              missing_text = "") %>%
  fmt_number(columns = 2:last_col(), 
             decimals = 5) 

gtsave(data = mosquito_dis_matrix , 
       filename = "plots/mosquito_dis_matrix.png")


#### preparing data for the Popart figure ----

nigripis_popart <- KBIMP_updatedspecies %>%
  filter(Species == "Aedes nigripes/impiger") %>%
  select(Sample, Sequence) 

nigripis_popart2 <- DNAStringSet(nigripis_popart$Sequence)

names(nigripis_popart2) <- nigripis_popart$Sample

alighned_nigripis_popart <- DNAStringSet(
  muscle::muscle(nigripis_popart2))

writeXStringSet(alighned_nigripis_popart,
                filepath = "processed-data/alighned_hexodontus.fasta",
                format = "fasta")

nigripis_popart_phydat <- as.phyDat(alighned_nigripis_popart, type = "DNA")

# Convert to character matrix
alighned_nigripis_popart_mat <- as.character(alighned_nigripis_popart)

# Write to NEXUS format (seqinr writes sequential by default)
write.nexus.data(alighned_nigripis_popart_mat, 
                 file = "processed-data/output_alignment_nigripis.nex"
                 ,format = "dna")







