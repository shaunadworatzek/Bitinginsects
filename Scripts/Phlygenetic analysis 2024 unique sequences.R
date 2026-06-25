library(Biostrings)


kbimp_mos_DNA_df <- KBIMP2024 %>%
  filter(Family %in% c("Culicidae", "Outgroup")) %>%
  select(Sample, Sequence) 
  
#preforming the alignment and converting to a PhyDat
  
kbimp_mos_DNA <- DNAStringSet(kbimp_mos_DNA_df$Sequence)
  
names(kbimp_mos_DNA) <- kbimp_mos_DNA_df$Sample
  
alighned_kbimpmos_DNA <- DNAStringSet(muscle::muscle(kbimp_mos_DNA))

browseSeqs(aligned_kbimpmos_DNA)

unique_seqsmos <- unique(alighned_kbimpmos_DNA)

writeXStringSet(unique_seqsmos,
                filepath = "processed-data/uniquemosseqforbold.fasta",
                format = "fasta")

kbimp_mos_phydat <- as.phyDat(unique_seqsmos, type = "DNA")
class(kbimp_mos_phydat) # is a "phyDat" object
length(kbimp_mos_phydat) # 68 has the  species as seen before 

BOLDID_mos <- read.csv(file = "processed-data/KBIMP2024_BOLDresults.csv")

BOLDID_mos2 <- BOLDID_mos  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[BOLD:(.*?)\\]",
          remove = FALSE) %>% 

  filter(!is.na(BOLDID)) %>%
  
  group_by(Query.ID) %>%
  
  slice_max(ID., n = 1, with_ties = FALSE) %>%
  
  ungroup() 

###### Building the tree ######

dist.mos <- dist.ml(kbimp_mos_phydat, ratio = TRUE, model = "JC69") 

#neighbor joining method

NJtree.kbimp.mos <- NJ(dist.mos)

plot(NJtree.kbimp.mos)

length(NJtree.kbimp.mos$tip.label) #68 tips as expected 

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

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.mos$node.label)

# Internal node IDs

internal_nodes <- (Ntip(rooted.bstree.mos)+1):(Ntip(rooted.bstree.mos)+Nnode(rooted.bstree.mos))


bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  mutate(bootstrap = (bootstrap *100)) %>%
  filter(bootstrap >= 60) 

##### Produce finalized tree #####

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

#produce the tree

mosquitotree2024 <- (ggtree(rooted.bstree.mos, layout = "rectangular", branch.length = TRUE) +
                       
                       geom_text(aes(label = BOLDID_mos2$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
                       
                       theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
                             legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, color = "red", fill = "white") 

mosquitotree2024

ggsave("plots/mosallseqtree.png", plot = mosquitotree2024, width = 5, height = 8, dpi = 300)


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

alighned_kbimpbf_DNA <- DNAStringSet(muscle::muscle(kbimp_bf_DNA))

BrowseSeqs(alighned_kbimpbf_DNA)

unique_seqsbf <- unique(alighned_kbimpbf_DNA)

writeXStringSet(unique_seqsbf,
                filepath = "processed-data/unique_bfseq_forbold.fasta",
                format = "fasta")

kbimp_bf_phydat <- as.phyDat(unique_seqsbf, type = "DNA")
class(kbimp_bf_phydat) # is a "phyDat" object
length(kbimp_bf_phydat) # 83 unique seq

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


BOLDID_bf <- read.csv(file = "processed-data/KBIMP2024_bfBOLD.csv")

BOLDID_bf2 <- BOLDID_bf  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[BOLD:(.*?)\\]",
          remove = FALSE) %>% 
  
  filter(!is.na(BOLDID)) %>%
  
  group_by(Query.ID) %>%
  
  slice_max(ID., n = 1, with_ties = FALSE) %>%
  
  ungroup() 

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.bf$node.label)

# Internal node IDs

internal_nodes <- (Ntip(rooted.bstree.bf)+1):(Ntip(rooted.bstree.bf)+Nnode(rooted.bstree.bf))


bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  mutate(bootstrap = (bootstrap *100)) %>%
  filter(bootstrap >= 60) 

##### Produce finalized tree #####

# trait data and the tree have the same order of species 

tree_tips_bf <- rooted.bstree.bf$tip.label

sample_name_bf <- BOLDID_bf2$Query.ID 

# checking if all species names are present in the tree

all(tree_tips_bf %in% sample_name_bf)

all(sample_name_bf %in% tree_tips_bf) 

setdiff(tree_tips_bf, sample_name_bf)

setdiff(sample_name_bf, tree_tips_bf)

#getting trait data set up as tip labels 

BOLDID_bf2 <- BOLDID_bf2[match(tree_tips_bf, BOLDID_bf2$Query.ID), ] 

BOLDID_bf2$tip_label_new <- paste0(BOLDID_bf2$BOLDID," (", BOLDID_bf2$ID., ")")


#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.bf)

#produce the tree

bftreespecies2024 <- (ggtree(rooted.bstree.bf, layout = "rectangular", branch.length = TRUE) +
                       
                       geom_text(aes(label = BOLDID_bf2$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
                       
                       theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
                             legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, color = "red", fill = "white") 

bftree2024

ggsave("plots/bftreespecies.png", plot = bftreespecies2024, width = 7, height = 10, dpi = 300)


ggsave("plots/bftreeboldid.png", plot = bftreespecies2024, width = 7, height = 10, dpi = 300)


##### species assignments for bf -----

BOLDIDspecies <- read_csv(file = "processed-data/BOLDIDspecies.csv")

# Combine BOLD tables
BOLDresults <- bind_rows(BOLDID_mos2, BOLDID_bf2)

KBIMP2024_updatedspecies2 <- KBIMP2024 %>%
  
  filter(Family %in% c("Culicidae", "Simuliidae")) %>%
  
  mutate(Species = str_replace(Species, "Aedes punctor", "Aedes punctor/Aedes hexodontus")) %>%
 
  left_join(BOLDresults, join_by("Sample" =="Query.ID")) %>%
  
  group_by(Sequence) %>%
  
 
  mutate(
    BOLDID = ifelse(
      any(!is.na(BOLDID)),        
      na.omit(BOLDID)[1], NA)) %>%
  
  ungroup() %>%
  
  left_join(BOLDIDspecies) %>%
  
  
  mutate(
    Species = ifelse(Species.x == "unknown" | is.na(Species.x), Species_BOLDID, Species.x),
    
    update_flag = case_when(
      Species.x == "unknown" | is.na(Species.x) ~ "Sequence Simularity",
      TRUE ~ "Probabolistic")) 
  


KBIMP2024_updatedspecies <- KBIMP2024_updatedspecies %>%
  filter(Family %in% c("Culicidae", "Simuliidae")) %>%
  select(Sample, Species)

KBIMP2024_updatedspecies2 <- KBIMP2024_updatedspecies2 %>%
  select(Sample, Species)


comparison <- KBIMP2024_updatedspecies %>%
  inner_join(KBIMP2024_updatedspecies2, by = "Sample", suffix = c("_df1", "_df2"))


differences <- comparison %>%
  filter(Species_df1 != Species_df2)

KBIMP2025 <- KBIMP2025 %>%
  mutate(Sample = str_replace(Sample, "^([A-Za-z]+)_(\\d+)", "\\1\\2")) %>%
  group_by(Sample) %>%
  mutate(Sample = paste0(Sample, "_", LETTERS[row_number()])) %>%
  ungroup() 


KBIMP <- bind_rows(KBIMP2024, KBIMP2025, Outgroup)

#### tree both years mosquitoes ---- 

kbimp_mos_DNA_df <- KBIMP %>%
  filter(Family %in% c("Culicidae", "Outgroup")) %>%
  select(Sample, Sequence) 

#preforming the alignment and converting to a PhyDat

kbimp_mos_DNA <- DNAStringSet(kbimp_mos_DNA_df$Sequence)

names(kbimp_mos_DNA) <- kbimp_mos_DNA_df$Sample

alighned_kbimpmos_DNA <- DNAStringSet(muscle::muscle(kbimp_mos_DNA))

browseSeqs(aligned_kbimpmos_DNA)

unique_seqsmos <- unique(alighned_kbimpmos_DNA)

writeXStringSet(unique_seqsmos,
                filepath = "processed-data/uniquemosseqforbold.fasta",
                format = "fasta")

kbimp_mos_phydat <- as.phyDat(unique_seqsmos, type = "DNA")
class(kbimp_mos_phydat) # is a "phyDat" object
length(kbimp_mos_phydat) # 68 has the  species as seen before

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

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.mos$node.label)

# Internal node IDs

internal_nodes <- (Ntip(rooted.bstree.mos)+1):(Ntip(rooted.bstree.mos)+Nnode(rooted.bstree.mos))


bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  mutate(bootstrap = (bootstrap *100)) %>%
  filter(bootstrap >= 60) 

##### Produce finalized tree #####

BOLDID_mos <- read.csv(file = "processed-data/BOLDID_aedes.csv")

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

#produce the tree

mosquitotree2024 <- (ggtree(rooted.bstree.mos, layout = "rectangular", branch.length = TRUE) +
                       
                       geom_text(aes(label = BOLDID_mos2$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
                       
                       theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
                             legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, color = "red", fill = "white") 

mosquitotree2024

ggsave("plots/mosallseqtree.png", plot = mosquitotree2024, width = 9, height = 16, dpi = 300)


#### both years black flies just Simulium genus ----

#filtering for Simuliidae, and defining the root as a chirmoidae sequence in the data set 

kbimp_bf_DNA_df <- KBIMP %>%
  filter(Genus %in% c("Simulium", "Outgroup")) %>%
  select(Sample, Sequence)  

# Convert the 'sequence' column to a DNAStringSet

kbimp_bf_DNA <- DNAStringSet(kbimp_bf_DNA_df$Sequence)

names(kbimp_bf_DNA) <- kbimp_bf_DNA_df$Sample

alighned_kbimpbf_DNA <- DNAStringSet(muscle::muscle(kbimp_bf_DNA))

BrowseSeqs(alighned_kbimpbf_DNA)

unique_seqsbf <- unique(alighned_kbimpbf_DNA)

writeXStringSet(unique_seqsbf,
                filepath = "processed-data/unique_bfseq_forbold.fasta",
                format = "fasta")

kbimp_bf_phydat <- as.phyDat(unique_seqsbf, type = "DNA")
class(kbimp_bf_phydat) # is a "phyDat" object
length(kbimp_bf_phydat) # 83 unique seq

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


BOLDID_bf <- read.csv(file = "processed-data/KBIMP_bfBOLD.csv")

BOLDID_bf2 <- BOLDID_bf  %>%
  
  extract(PID..BIN.,
          into = c("ID", "BOLDID"),
          regex = "(.*?)\\[BOLD:(.*?)\\]",
          remove = FALSE) %>% 
  
  filter(!is.na(BOLDID)) %>%
  
  group_by(Query.ID) %>%
  
  slice_max(ID., n = 1, with_ties = FALSE) %>%
  
  ungroup() 

# Numeric bootstrap values for plotting on tree

bs_numeric <- as.numeric(rooted.bstree.bf$node.label)

# Internal node IDs

internal_nodes <- (Ntip(rooted.bstree.bf)+1):(Ntip(rooted.bstree.bf)+Nnode(rooted.bstree.bf))


bs_tibble <- tibble(node = internal_nodes, bootstrap = bs_numeric) %>%
  mutate(bootstrap = (bootstrap *100)) %>%
  filter(bootstrap >= 60) 

##### Produce finalized tree #####

# trait data and the tree have the same order of species 

tree_tips_bf <- rooted.bstree.bf$tip.label

sample_name_bf <- BOLDID_bf2$Query.ID 

# checking if all species names are present in the tree

all(tree_tips_bf %in% sample_name_bf)

all(sample_name_bf %in% tree_tips_bf) 

setdiff(tree_tips_bf, sample_name_bf)

setdiff(sample_name_bf, tree_tips_bf)

#getting trait data set up as tip labels 

BOLDID_bf2 <- BOLDID_bf2[match(tree_tips_bf, BOLDID_bf2$Query.ID), ] 

BOLDID_bf2$tip_label_new <- paste0(BOLDID_bf2$BOLDID," (", BOLDID_bf2$ID., ")")


#plotting the tree with trait data 

node <- 1:Ntip(rooted.bstree.bf)

#produce the tree

bftreespecies2024 <- (ggtree(rooted.bstree.bf, layout = "rectangular", branch.length = TRUE) +
                        
                        geom_text(aes(label = BOLDID_bf2$tip_label_new[node]), hjust = -0.05, size =4, fontface = "italic") +
                        
                        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
                              legend.position = "top")) %<+% 
  
  bs_tibble +
  geom_label2(aes(label = bootstrap), hjust = 0.7, size = 3, color = "red", fill = "white") 

bftree2024

ggsave("plots/bftreespecies2.png", plot = bftreespecies2024, width = 9, height = 17, dpi = 300)


ggsave("plots/bftreeboldid.png", plot = bftreespecies2024, width = 7, height = 10, dpi = 300)

