
KBIMP2024 <- read_tsv(file = "processed-data/KBIMP2024_filteredCOI.tsv")
KBIMP2025 <- read_tsv(file = "processed-data/KBIMP2025_filteredCOI.tsv")
Outgroup <- read_csv(file = "raw-data2/Outgroup.csv")

KBIMP2024 <- KBIMP2024 %>%
  full_join(kbimp2024_sampledata, join_by(Sample == SampleID)) %>%
  
  mutate(Sample = paste0(FieldID,case_when(
                             Family == "Simuliidae" ~ "_BF",
                             Family == "Culicidae" ~ "_M", TRUE ~ ""))) %>%
  select(Sample, Read_Count, Kingdom, Phylum, Class, Order, Family, Genus, Species, Probability_Genus, Probability_Species, Sequence) %>%
  group_by(Sample) %>%
  mutate(Sample = paste0(Sample, "_", LETTERS[row_number()])) %>%
  ungroup()


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

###### Building the tree ######

dist.mos <- dist.ml(kbimp_mos_phydat, ratio = TRUE, model = "JC69") 

#neighbor joining method

NJtree.kbimp.mos <- NJ(dist.mos)

plot(NJtree.kbimp.mos)

length(NJtree.kbimp.mos$tip.label) #68 tips as expected 

# Fit the initial tree using a simple pml 

pml.tree.kbimp.mos <- pml(NJtree.kbimp.mos, kbimp_mos_phydat, k = 4, model = "JC", method = "unrooted")

plot(pml.tree.kbimp.mos$tree) 

pml.tree.kbimp.mosop <- pml_bb(pml.tree.kbimp.mos, model = "JC")

#bootstrapping to determine reliability

bs <- bootstrap.pml(pml.tree.kbimp.mosop, bs = 1000, optNni = TRUE, multicore = TRUE)

tree_with_bs <- plotBS(pml.tree.kbimp.mos$tree, bs)

#rooting the tree

rooted.bstree.mos <- root(tree_with_bs, outgroup = "Outgroup", resolve.root = TRUE)

plot(rooted.bstree.mos)

tree_with_bs <- plotBS(rooted.bstree.mos, bs)

##### preparing labels for finalzed tree #####

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.mos$node.label)

internal_nodes <- (Ntip(rooted.bstree.mos)+1):(Ntip(rooted.bstree.mos)+Nnode(rooted.bstree.mos))

bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  
  mutate(bootstrap = (bootstrap *100)) %>%
  
  filter(bootstrap >= 60) 

# brining in BOLD data for tree 

BOLDID_mos <- read.csv(file = "processed-data/BOLDaedescombined.csv")

BOLDID_mos2 <- BOLDID_mos  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[BOLD:(.*?)\\]",
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


#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.mos)

##### Produce finalized tree #####

mosquitotree2024 <- (ggtree(rooted.bstree.mos, layout = "rectangular", branch.length = TRUE) +
                       
                       geom_text(aes(label = BOLDID_mos2$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
                       
                       theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
                             legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, color = "red", fill = "white") 

mosquitotree2024

ggsave("plots/mosallseqtree.png", plot = mosquitotree2024, width = 9, height = 16, dpi = 300)

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

rooted.bstree.sim <- root(tree_with_bssim, outgroup = "Outgroup", resolve.root = TRUE)

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
          regex = "(.*?)\\[BOLD:(.*?)\\]",
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
                        
                        geom_text(aes(label = BOLDID_sim2$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
                        
                        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
                              legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, color = "red", fill = "white") 

simtreespecies2024 

ggsave("plots/Simulium_treespecies.png", plot = simtreespecies2024, width = 15, height = 17, dpi = 300)



#### both years black flies not Simulium genus ----

kbimp_bf_DNA_df <- KBIMP %>%
  filter(Family %in% c("Simuliidae", "Outgroup")) %>%
  filter(!Genus == "Simulium") %>%
  select(Sample, Sequence)  

##### aligning and preparing phydat #####

kbimp_bf_DNA <- DNAStringSet(kbimp_bf_DNA_df$Sequence)

names(kbimp_bf_DNA) <- kbimp_bf_DNA_df$Sample

alighned_kbimpbf_DNA <- DNAStringSet(muscle::muscle(kbimp_bf_DNA))

#BrowseSeqs(alighned_kbimpbf_DNA)

unique_seqsbf <- unique(alighned_kbimpbf_DNA)

writeXStringSet(unique_seqsbf,
                filepath = "processed-data/unique_bfseq_forbold.fasta",
                format = "fasta")

kbimp_bf_phydat <- as.phyDat(unique_seqsbf, type = "DNA")
class(kbimp_bf_phydat) # is a "phyDat" object
length(kbimp_bf_phydat) # 83 unique seq

###### Building the tree ######

#create a new dist matrix

dist.medoid.bf <- dist.ml(kbimp_bf_phydat, ratio = TRUE, model = "JC69") 

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

##### preparing labels for finalzed tree #####

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.bf$node.label)

internal_nodes <- (Ntip(rooted.bstree.bf)+1):(Ntip(rooted.bstree.bf)+Nnode(rooted.bstree.bf))

bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  mutate(bootstrap = (bootstrap *100)) %>%
  filter(bootstrap >= 60) 

#bringing in bold data for final tree 

BOLDID_bfnotsim <- read.csv(file = "processed-data/BOLDID_notsimulidae.csv")

BOLDID_bfnotsim2 <- BOLDID_bfnotsim  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[BOLD:(.*?)\\]",
          remove = FALSE) %>% 
  
  filter(!is.na(BOLDID)) %>%
  
  filter(ID. >= 94) %>%
  
  group_by(Query.ID) %>%
  
  slice_max(ID., n = 1, with_ties = FALSE) %>%
  
  ungroup() 

# trait data and the tree have the same order of species 

tree_tips_bf <- rooted.bstree.bf$tip.label

sample_name_bf <- BOLDID_bfnotsim2$Query.ID 

# checking if all species names are present in the tree

all(tree_tips_bf %in% sample_name_bf)

all(sample_name_bf %in% tree_tips_bf) 

setdiff(tree_tips_bf, sample_name_bf)

setdiff(sample_name_bf, tree_tips_bf)

#getting trait data set up as tip labels 

BOLDID_bfnotsim2 <- BOLDID_bfnotsim2[match(tree_tips_bf, BOLDID_bfnotsim2$Query.ID), ] 

BOLDID_bfnotsim2$tip_label_new <- paste0(BOLDID_bfnotsim2$BOLDID," (", BOLDID_bfnotsim2$ID., ")")


#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.bf)

##### Produce finalized tree #####

bftreespecies2024 <- (ggtree(rooted.bstree.bf, layout = "rectangular", branch.length = TRUE) +
                        
                        geom_text(aes(label = BOLDID_bfnotsim2$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
                        
                        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
                              legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, color = "red", fill = "white") 

bftreespecies2024 

ggsave("plots/NotSimulium_treespecies.png", plot = bftreespecies2024, width = 28, height = 17, dpi = 300)


##### species assignments and comparing to mediod method -----

BOLDIDspecies <- read_csv(file = "processed-data/BOLDIDspecies.csv")

# Combine BOLD tables
BOLDresults <- bind_rows(BOLDID_mos2, BOLDID_bf2, BOLDID_bfnotsim2)

BOLDresults <- BOLDresults %>%
  filter(!Query.ID == "Outgroup") %>%
  filter(!is.na(Query.ID)) %>%
  select(Query.ID, BOLDID)

KBIMP_updatedspecies2 <- KBIMP %>%
  filter(Family %in% c("Culicidae", "Simuliidae")) %>%
  mutate(Species = str_replace(Species, "Aedes punctor", "Aedes punctor/Aedes hexodontus")) %>%
  left_join(BOLDresults, join_by("Sample" =="Query.ID")) %>%
  group_by(Sequence) %>%
  mutate(BOLDID = ifelse(any(!is.na(BOLDID)), na.omit(BOLDID)[1], NA)) %>%
  ungroup() %>%
  left_join(BOLDIDspecies) %>%
  mutate(update_flag = case_when(
      Species == "unknown" | is.na(Species) ~ "Sequence Simularity",
      TRUE ~ "Probabolistic")) %>%
    mutate(Species = coalesce(Species_BOLDID, Species)) %>%
    select(-Species_BOLDID) %>%
    filter(!Species == "unknown") %>%
  select(Sample, Species, Family, update_flag) %>%
  mutate(Sample = str_remove(Sample, "_.*")) %>%
  filter(!is.na(Species)) %>% 
  distinct(across(-update_flag), .keep_all = TRUE)

write_tsv(KBIMP_updatedspecies2, "processed-data/KBIMP_updatedspecies.tsv")
