
#loading required packages 

library(tidyr)
library(readr)
library(viridis)
library(Biostrings) #for DNA alignment and further analysis 
library(ape)      #for centroid determination 
library(muscle) #for sequence alignment
library(phangorn) #for tree construction 
library(parallel) #for running bootstrapping in parallel
library(ggplot2)    #for plotting 
library(ggtree)     #for plotting phylogeny trees 
library(seqinr)
library(DECIPHER)

#setting working directory to coursefile location
setwd("C:/Users/sdwor/OneDrive - University of Guelph/1- MASTERS DEGREE/Bitinginsects/Source_files")

#### PART 1 - Mosquito phylogeny ----

#opening fasta file with just the mosquito sequences 

CBAY_mos_DNA <- readDNAStringSet(file = "../Data/CBAY2024mosalignment.fasta")
class(CBAY_mos_DNA)

#preforming an alignment 

alighned_mos_DNA <- DNAStringSet(muscle::muscle(CBAY_mos_DNA))

#converting fasta file into a phyDat format 

CBAY_mos_phydat <- as.phyDat(alighned_mos_DNA, type = "DNA")
class(CBAY_mos_phydat) # is a "phyDat" object
length(CBAY_mos_phydat) # has the 1272 species as seen before

#determining the best model for the ML tree 

#modelTest <- modelTest(CBAY_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) # JC is the best model because it has the lowest AIC, second lowest BIC and high loglik


# making a tree with maximum likelihood using the JC model
# creating a distance matrix
dist.CBAY.mos <- dist.ml(CBAY_mos_phydat, ratio = TRUE, model = "JC69") 

# creating a tree using the neighbor joining method
NJtree.CBAY.mos <- NJ(dist.CBAY.mos)
plot(NJtree.CBAY.mos) #a tree was produced
length(NJtree.CBAY.mos$tip.label) #1272 tips as expected 

# Fit the initial tree using a simple pml 
pml.tree.CBAY.mos <- pml(NJtree.CBAY.mos, CBAY_mos_phydat, k = 4, model = "GTR+I", method = "unrooted")
plot(pml.tree.CBAY.mos$tree) #a tree was produced
pml.tree.CBAY.mosop <- pml_bb(pml.tree.CBAY.mos, model = "GTR+I")

bs <- bootstrap.pml(pml.tree.CBAY.mosop, bs = 1000, optNni = TRUE, multicore = TRUE)
tree_with_bs <- plotBS(pml.tree.CBAY.mos$tree, bs)

#rooting the tree
rooted.bstree.mos <- root(tree_with_bs, outgroup = "albopitcus_root", resolve.root = TRUE)
plot(rooted.bstree.mos)

#adding in a label for the root of the tree 

newrow <- data.frame(SampleID ="albopitcus_root", ExactSite = "Root")

#comb ing aligned sequences with the metadata for sites (site metadata is made in ecological analysis code)

CBAY_mosmetadata <- CBAY2024_sites %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, SamplingProtocol, FieldID, ExactSite, Lat, Lon, Genus) %>%
  filter(Genus == "Aedes") %>%
  bind_rows(newrow)
  

#making a tree with the sampling locations labels on the tree tips 

# the next few steps are done to make sure that the trait data and the tree have the same order of species for the tree construction
tree_tips <- rooted.bstree.mos$tip.label 
sample_name <- CBAY_mosmetadata$SampleID # getting the sample names from the trait data
# checking if all species names are present in the tree
all(tree_tips %in% sample_name) #FALSE not all sample_names in tree tips
all(sample_name %in% tree_tips) #TRUE all the sample names are in the three tips
setdiff(tree_tips, sample_name)
setdiff(sample_name, tree_tips)
CBAY_mosmetadata$SampleID <- CBAY_mosmetadata[match(tree_tips, sample_name), ] # reordering the trait data to match the order of the tree tip labels

#preparing trait data for phylo community plots 
length(unique(CBAY_mosmetadata$ExactSite)) 

# Assign colors to groups
group_colors <- c("dew line rd"= "#00D400", "freshwater river" = "#FF2A7F",
                  "Grenier lake"= "#008066", "long pt creek" ="#FF6600", "gravel pit" ="#0068FF")
# Add colors to the trait data
CBAY_mosmetadata$color <- group_colors[CBAY_mosmetadata$ExactSite]

#plotting the tree with trait data 
node <- 1:Ntip(rooted.bstree.mos)



nicetree <- (
  ggtree(rooted.bstree.mos, layout = "dendrogram", branch.length = TRUE) +
    geom_tippoint(aes(color = CBAY_mosmetadata$ExactSite[node]), size = 2) +
    scale_color_manual(values = group_colors, name = "Sampling Location") +
    theme_tree2(legend.position = "right") +
    theme(
      legend.title = element_text(size = 12),      
      legend.text = element_text(size = 10)
    )
)
nicetree
ggsave("mosphylo.png", plot = nicetree, width = 6, height = 6, dpi = 300)

#### - PART 2 - making the data file for the haplotype network  ----

# Hierarchical clustering
hc <- hclust(dist.CBAY.mos)
groups <- cutree(hc, k = 3)

# View assignments
table(groups)

seqs_group1 <- alighned_mos_DNA[names(groups)[groups == 1]]
seqs_group2 <- alighned_mos_DNA[names(groups)[groups == 2]]

names(seqs_group1) <- gsub("_", "", names(seqs_group1))
names(seqs_group1) <- gsub("-", "", names(seqs_group1))



CBAY_mosmetadata$SampleID <- gsub("_", "", CBAY_mosmetadata$SampleID)
CBAY_mosmetadata$SampleID <- gsub("-", "", CBAY_mosmetadata$SampleID)

CBAY_mosmetadatagroup1 <- CBAY_mosmetadata %>%
  filter(SampleID %in% names(seqs_group1)) %>%
  select(SampleID, ExactSite) %>%
  column_to_rownames(var = "SampleID")


# Convert to binary matrix
binary_matrix_traitg1 <- model.matrix(~ ExactSite - 1, data = CBAY_mosmetadatagroup1)

binary_matrix_traitg1 <- binary_matrix_traitg1 %>% as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  select("SampleID","ExactSitedew line rd" , "ExactSitefreshwater river","ExactSiteGrenier lake",
         "ExactSitelong pt creek", "ExactSitegravel pit")

write_tsv(binary_matrix_traitg1, "mosmetadatag1.tsv")


# 2. Compare exact entries side by side
head(rownames(seqs_group1), 10)
head(CBAY_mosmetadatagroup1$SampleID, 10)

# 3. Check for whitespace
any(grepl("^\\s|\\s$", CBAY_mosmetadatagroup2$SampleID))  # leading/trailing in trait
any(grepl("^\\s|\\s$", rownames(seqs_group2)))            # leading/trailing in seqs

# 4. Remove all whitespace and compare again
sampleIDs_clean <- trimws(CBAY_mosmetadatagroup2$SampleID)
seqnames_clean <- trimws(rownames(seqs_group2))

all(sampleIDs_clean %in% seqnames_clean)
all(seqnames_clean %in% sampleIDs_clean)

# 5. See the actual mismatches
setdiff(sampleIDs_clean, seqnames_clean)
setdiff(seqnames_clean, sampleIDs_clean)

CBAY_mosmetadatagroup2$SampleID <- trimws(CBAY_mosmetadatagroup2$SampleID)
rownames(seqs_group2) <- trimws(rownames(seqs_group2))

#making a nexus file 

write.nexus.data(seqs_group1, file = "output_alignment.nex")


all(rownames(CBAY_mosmetadatagroup1) %in% seqs_group1$ranges$NAMES)
all(seqs_group1$ranges$NAMES %in% rownames(CBAY_mosmetadatagroup1))
setdiff(rownames(CBAY_mosmetadatagroup1), seqs_group1$ranges$NAMES)
setdiff(seqs_group1$ranges$NAMES, CBAY_mosmetadatagroup1$SampleID)


names(seqs_group2) <- gsub("_", "", names(seqs_group2))
names(seqs_group2) <- gsub("-", "", names(seqs_group2))

CBAY_mosmetadatagroup2 <- CBAY_mosmetadata %>%
  filter(SampleID %in% names(seqs_group2)) %>%
  select(SampleID, ExactSite) %>%
  column_to_rownames(var = "SampleID")

all(rownames(CBAY_mosmetadatagroup2) %in% seqs_group2$ranges$NAMES)
all(seqs_group2$ranges$NAMES %in% rownames(CBAY_mosmetadatagroup2))
setdiff(rownames(CBAY_mosmetadatagroup2), seqs_group2$ranges$NAMES)
setdiff(seqs_group2$ranges$NAMES, CBAY_mosmetadatagroup2$SampleID)

# Convert to binary matrix
binary_matrix_traitg2 <- model.matrix(~ ExactSite - 1, data = CBAY_mosmetadatagroup2)

binary_matrix_traitg2 <- binary_matrix_traitg2 %>% as.data.frame() %>%
  rownames_to_column(var = "SampleID")


write_tsv(binary_matrix_traitg2, "mosmetadatag2.tsv")

#making a nexus file 

write.nexus.data(seqs_group2, file = "output_alignment2.nex")


#### PART 3 - looking at all of the black fly samples ---- 

CBAY_bf_DNAall <- readDNAStringSet(file = "../Data/CBAY_bf_2024.fasta")
class(CBAY_bf_DNAall)

#preforming an alignment 

alighned_bf_DNA <- DNAStringSet(muscle::muscle(CBAY_bf_DNAall))

#BrowseSeqs(alighned_bf_DNA)

#converting fasta file into a phyDat format 

CBAY_bf_phydat <- as.phyDat(alighned_bf_DNA, type = "DNA")
class(CBAY_bf_phydat) # is a "phyDat" object
length(CBAY_bf_phydat)

#determining the best model for the ML tree 

#modelTest <- modelTest(CBAY_bf_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) # JC is the best model because it has the lowest AIC, second lowest BIC and high loglik

# making a tree with maximum likelihood using the JC model
# creating a distance matrix
dist.CBAY.bf <- dist.ml(CBAY_bf_phydat, ratio = TRUE, model = "JC69") 

# creating a tree using the neighbor joining method
NJtree.CBAY.bf <- NJ(dist.CBAY.bf)
plot(NJtree.CBAY.bf, cex =0.5) #a tree was produced
length(NJtree.CBAY.bf$tip.label) #1272 tips as expected 

# Fit the initial tree using a simple pml 
pml.tree.CBAY.bf <- pml(NJtree.CBAY.bf, CBAY_bf_phydat, k = 4, model = "GTR+I", method = "unrooted")
plot(pml.tree.CBAY.bf$tree) #a tree was produced
pml.tree.CBAY.bfop <- pml_bb(pml.tree.CBAY.bf, model = "GTR+I")

bsbf <- bootstrap.pml(pml.tree.CBAY.bfop, bs = 1000, optNni = TRUE, multicore = TRUE)
tree_with_bs <- plotBS(pml.tree.CBAY.bf$tree, bsbf)  # Only show support ≥ 50%

rooted.bstree.bf <- root(tree_with_bs, outgroup = "root\t", resolve.root = TRUE, keep.root.edge = TRUE)
is.rooted(rooted.bstree.bf)
plot(rooted.bstree.bf)

#adding a label for the root 

newrowbf <- data.frame(SampleID ="root\t", Species = "Root")

#combing with metadata 

CBAY_bfmetadata <- CBAY2024_sites %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  filter(Family == "Simuliidae") %>% 
  select(SampleID, Species)%>%
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum*", Species)) %>%
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare*", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium sp. Unclassified", Species)) %>%
  filter(Species != "unknown") %>%
  bind_rows(newrowbf)

#making a tree with the sampling locations labels on the tree tips 

# the next few steps are done to make sure that the trait data and the tree have the same order of species for the tree construction
tree_tips <- rooted.bstree.bf$tip.label # getting the tip labels from the tree
sample_name <- CBAY_bfmetadata$SampleID # getting the sample names from the trait data
# checking if all species names are present in the tree
all(tree_tips %in% sample_name) #FALSE not all sample_names in tree tips
all(sample_name %in% tree_tips) #TRUE all the sample names are in the three tips 
CBAY_bfmetadata <- CBAY_bfmetadata[match(tree_tips, CBAY_bfmetadata$SampleID), ] # reordering the trait data to match the order of the tree tip labels
setdiff(tree_tips, sample_name)
setdiff(sample_name, tree_tips)
#preparing trait data for phylo community plots 
length(unique(CBAY_bfmetadata$Species)) #13 orders present in the dataset



# Assign colors to groups
group_colors <- setNames(c("#D62728",  
                           "#F0E442",
                           "#FF5733",  
                           "#CC79A7",
                           "#E69F00",  
                           "#0072B2",  
                           "#009E73",
                           "#FF69B4", "#8A2BE2", "#7AC5CD", "grey45"), unique(CBAY_bfmetadata$Species))
# Add colors to the trait data
CBAY_bfmetadata$color <- group_colors[CBAY_bfmetadata$Species]

#plotting the tree with trait data 
node <- 1:Ntip(rooted.bstree.bf)


nicetree <- (
  ggtree(rooted.bstree.bf, layout = "dendrogram", branch.length = TRUE, root.position = 0) +
    geom_tippoint(aes(color = CBAY_bfmetadata$Species[node]), size = 1.3) +
    scale_color_manual(values = group_colors, name = "Species") +
    theme_tree2(legend.position = "left") +
    theme(
      legend.title = element_text(size = 10, face = "bold"),      
      legend.text = element_text(size = 8)
    )
)
plot(nicetree) 
ggsave("phylogeny.png", plot = nicetree, width = 6, height = 3, dpi = 300)

#### - trying to make a haplotype nwtrok for the subpsillium  ----

sampleidsupbsillium <- CBAY_bfmetadata  %>%
  filter(Species %in% c("Simulium sp. Unclassified", "Simulium subpusillum"))

fasta_supb <- alighned_bf_DNA[names(alighned_bf_DNA) %in% sampleidsupbsillium$SampleID]

names(fasta_supb) <- gsub("_", "", names(fasta_supb))
names(fasta_supb) <- gsub("-", "", names(fasta_supb))

metabfdatasim$SampleID <- gsub("_", "", metabfdatasim$SampleID)
metabfdatasim$SampleID <- gsub("-", "", metabfdatasim$SampleID)

CBAY_bfmetasupb <- metabfdatasim %>%
  filter(SampleID %in% names(fasta_supb)) %>%
  select(SampleID, ExactSite) %>%
  column_to_rownames(var = "SampleID")


# Convert to binary matrix
binary_matrix_traitsupb <- model.matrix(~ ExactSite - 1, data = CBAY_bfmetasupb)

binary_matrix_traitsupb <- binary_matrix_traitsupb %>% as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  select("SampleID","ExactSitedew line rd" , "ExactSitefreshwater river","ExactSiteGrenier lake",
         "ExactSitelong pt creek")

write_tsv(binary_matrix_traitsupb, "binary_matrix_traitsupb.tsv")


# 2. Compare exact entries side by side
head(names(fasta_supb), 10)
head(binary_matrix_traitsupb$SampleID, 10)

# 3. Check for whitespace
any(grepl("^\\s|\\s$", binary_matrix_traitsupb$SampleID))  # leading/trailing in trait
any(grepl("^\\s|\\s$", names(fasta_supb)))            # leading/trailing in seqs

# 4. Remove all whitespace and compare again
sampleIDs_clean <- trimws(binary_matrix_traitsupb$SampleID)
seqnames_clean <- trimws(names(fasta_supb))

all(sampleIDs_clean %in% seqnames_clean)
all(seqnames_clean %in% sampleIDs_clean)

#making a nexus file 

write.nexus.data(fasta_supb, file = "fasta_supb.nex")

#### - trying to make a haplotype nwtrok for the simuluum group  ----

sampleidsimulium <- CBAY_bfmetadata  %>%
  filter(Species %in% c("Simulium noelleri", "Simulium decorum*", "Simulium vulgare*", "Simulium decimatum"))

fasta_simiulum <- alighned_bf_DNA[names(alighned_bf_DNA) %in% sampleidsupbsimulium$SampleID]

names(fasta_simiulum) <- gsub("_", "", names(fasta_simiulum))
names(fasta_simiulum) <- gsub("-", "", names(fasta_simiulum))

metabfdatasim$SampleID <- gsub("_", "", metabfdatasim$SampleID)
metabfdatasim$SampleID <- gsub("-", "", metabfdatasim$SampleID)

CBAY_bfmetasimulium <- metabfdatasim %>%
  filter(SampleID %in% names(fasta_simiulum)) %>%
  select(SampleID, ExactSite) %>%
  column_to_rownames(var = "SampleID")


# Convert to binary matrix
binary_matrix_traitsimulium <- model.matrix(~ ExactSite - 1, data = CBAY_bfmetasimulium)

binary_matrix_traitsimulium  <- binary_matrix_traitsimulium  %>% as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  select("SampleID","ExactSitedew line rd" , "ExactSitefreshwater river",
         "ExactSitelong pt creek")

write_tsv(binary_matrix_traitsimulium, "binary_matrix_traitsimuliilium.tsv")


# 2. Compare exact entries side by side
head(names(fasta_simiulum), 10)
head(binary_matrix_traitsimulium$SampleID, 10)

# 3. Check for whitespace
any(grepl("^\\s|\\s$", binary_matrix_traitsimulium$SampleID))  # leading/trailing in trait
any(grepl("^\\s|\\s$", names(fasta_simiulum)))            # leading/trailing in seqs

#making a nexus file 

write.nexus.data(fasta_simiulum, file = "fasta_simiulum.nex")



