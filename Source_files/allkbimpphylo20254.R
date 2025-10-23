
#loading required packages 


library(tidyverse)
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

#### - Preparing the site metadata for analysis ----

kbimp2024_sampledata <- read_csv(file = "../raw-data/KBIMP2024_specimendata.csv")

view(kbimp2024_sampledata)

KBIMP_2024_siteslook <- kbimp2024_sampledata %>%
  group_by(ExactSite) %>%
  summarise(countsites = n_distinct(FieldID)) 


ggplot(KBIMP_2024_siteslook, aes(y= countsites, x= ExactSite)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(KBIMP_2024_siteslook)

#fixing all the messed up site names

kbimp2024_sampledata <- kbimp2024_sampledata %>%
  mutate(ExactSite = tolower(ExactSite)) %>%
  mutate(ExactSite = gsub("by the creek at ", "", ExactSite)) %>%
  mutate(ExactSite = gsub("on ", "", ExactSite)) %>%
  mutate(ExactSite = gsub("4 mile trail", "trail to 4 mile", ExactSite)) %>%
  mutate(ExactSite = gsub("trail to 4 mile bay", "trail to 4 mile", ExactSite)) %>%
  mutate(ExactSite = gsub("4 mile road", "trail to 4 mile", ExactSite)) %>%
  mutate(ExactSite = gsub("by pond", "", ExactSite)) %>%
  mutate(ExactSite = gsub("road", "rd", tolower(ExactSite))) %>%
  mutate(ExactSite = gsub(" resample", "", ExactSite)) %>%
  mutate(ExactSite = gsub("- site 1", "", ExactSite)) %>%
  mutate(ExactSite = gsub("stream ", "", ExactSite)) %>%
  mutate(ExactSite= gsub("west arm/gravel pit near sea", "gravel pit", ExactSite)) %>%
  mutate(ExactSite= gsub("first", "1st", ExactSite)) %>%
  mutate(ExactSite = gsub("\\s+", " ", ExactSite)) %>%
  mutate(ExactSite= gsub("1st creek/", "", ExactSite)) %>%
  mutate(ExactSite= gsub("kitigak river", "Grenier lake", ExactSite)) %>%
  mutate(ExactSite= gsub("2nd culvert pelly rd", "Grenier lake", ExactSite)) %>%
  mutate(ExactSite= gsub("into grenier lake", "Grenier lake", ExactSite)) %>%
  mutate(ExactSite= gsub("between 1st lake and grenier", "Grenier lake", ExactSite)) %>%
  mutate(ExactSite= gsub("inland area", "freshwater river", ExactSite)) %>%
  mutate(ExactSite= gsub("creek south of pelly rd", "freshwater river", ExactSite)) %>%
  mutate(ExactSite= gsub("s of pelly rd near swimming hole", "freshwater river by picnic area", ExactSite)) %>%
  mutate(ExactSite = str_replace(ExactSite, "\\s+$", "")) %>%
  mutate(ExactSite= gsub("freshwater river by picnic area", "freshwater river", ExactSite)) %>%
  mutate(ExactSite= gsub("the ravine", "long pt creek", ExactSite)) %>%
  mutate(ExactSite= gsub("1st creek", "long pt creek", ExactSite)) %>%
  mutate(ExactSite= gsub("past old town", "freshwater river", ExactSite)) %>%
  mutate(ExactSite= gsub("dew line rd freshwater creek", "freshwater river", ExactSite)) %>%
  mutate(ExactSite= gsub("freshwater river pelly rd", "freshwater river", ExactSite)) %>%
  mutate(ExactSite= gsub("creek by long pt beach", "dew line rd", ExactSite)) %>%
  mutate(ExactSite = gsub("road", "trail", ExactSite)) %>%
  mutate(ExactSite = str_trim(ExactSite, side = "both"))

KBIMP_2024_siteslook <- kbimp2024_sampledata %>%
  group_by(ExactSite) %>%
  summarise(countsites = n_distinct(FieldID)) 


ggplot(KBIMP_2024_siteslook, aes(y= countsites, x= ExactSite)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(KBIMP_2024_siteslook)

#fixing the samples which were switches during extractions

kbimp2024_sampledata <- kbimp2024_sampledata %>%
  mutate(SampleID= gsub("KBIMP_004_H11", "KBIMP-_006_G3", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_D12", "KBIMP-_006_G4", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_E12", "KBIMP-_006_G5", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_F12", "KBIMP-_006_G6", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_G12", "KBIMP-_006_G7", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_H10", "KBIMP-_006_G2", SampleID)) 


#### - Filtering seq data for controls and combining different sequence runs ----

#input raw sequence data 

KGLKTK_COI_Domcontigs <- read_tsv(file = "../raw-data/Shauna_COI_KGLKTK2024_TaxonomicAssignments_DominantContigs.tsv")
CBAYCOIrestplates_Domcontigs <- read_tsv(file = "../raw-data/Shauna_COI_CBAY2024_TaxonomicAssignments_DominantContigs.tsv")
CBAYCOIplate3_Domcontigs <- read_tsv(file = "../raw-data/Shauna_100dilution_TaxonomicAssignments_DominantContigs.tsv")

View(KGLKTK_COI_Domcontigs)
View(CBAYCOIrestplates_Domcontigs)
View(CBAYCOIplate3_Domcontigs)

length(KGLKTK_COI_Domcontigs$Sample) #168
length(CBAYCOIrestplates_Domcontigs$Sample) #439
length(CBAYCOIplate3_Domcontigs$Sample) #100

#Filtering out/ accoutning for negative controls 
#negative control on plate 3

CBAYCOIplate3_Domcontigs [100, "Sample"] <- "Control"
 
CBAYCOIplate3_Domcontigs <- CBAYCOIplate3_Domcontigs %>%
  filter(!grepl("10dil", Sample)) %>%
  mutate(ReadCount = ReadCount - ReadCount[Sample == 'Control']) %>%
  filter(!ReadCount <= 0) 
  

#negative controls on plate 1,2,4,5,6

CBAYCOIrestplates_Domcontigs <- CBAYCOIrestplates_Domcontigs %>%
  filter(!Sample %in% c("KBIMP_004_H11", "KBIMP_004_D12", "KBIMP_004_E12", 
                        "KBIMP_004_F12", "KBIMP_004_G12", "KBIMP_004_H10", "KBIMP-_006_H1",
                        "KBIMP-_006_H10", "KBIMP-_006_H4", "KBIMP-_006_H5", "KBIMP-_006_H6",
                        "KBIMP-_006_H8")) %>%
  add_row(Plate ='KBIMP_01_CBAY2024', Sample = 'KBIMP_01_H12', ReadCount = 0) %>%
  group_by(Plate) %>%
  mutate(ReadCount = ReadCount - ReadCount[grepl("H12", Sample)]) %>%
  ungroup() %>%
  filter(!ReadCount <= 0)


#negative controls for Kugluktuk plate 

KGLKTK_COI_Domcontigs <- KGLKTK_COI_Domcontigs %>%
  mutate(ReadCount = ReadCount - ReadCount[Sample == 'KBIMP-_007_H12']) %>%
  filter(!ReadCount <= 0)

length(KGLKTK_COI_Domcontigs$Sample) #120 lost 48
length(CBAYCOIrestplates_Domcontigs$Sample) #354 lost 76 (but this includes the 7 repeats)
length(CBAYCOIplate3_Domcontigs$Sample) #95 lost 0 (there were doubles of some samples)

#combining plate 3 with the rest of the plates

CBAY_Domcontigs <- rbind(CBAYCOIrestplates_Domcontigs, CBAYCOIplate3_Domcontigs)

View(CBAY_Domcontigs)

KBIMP_Domcontigs <- rbind(CBAY_Domcontigs, KGLKTK_COI_Domcontigs)

write_tsv(KBIMP_Domcontigs, "../processed-data/alldata_2024_domcontigs.tsv")

#### - Investigating data ----

#lets look at the amount of wrong IDs (or IDs that aren't black flies/mosquitoes)

#looking at it by Family
KBIMP_Domcontigs_family <- KBIMP_Domcontigs %>%
  select(Family) %>%
  count(Family, name = "countfamily")

view(KBIMP_Domcontigs_family)

ggplot(KBIMP_Domcontigs_family, aes(y=countfamily, x=Family)) +
  geom_col(fill = "skyblue3") +
  theme_bw()  #7 families which weren't mosquitoes or black flies 

rm(KBIMP_Domcontigs_family)

#now lets look at the genus diversity of only black flies and mosquitoes

KBIMP_Domcontigs_bitinggenus <- KBIMP_Domcontigs %>%
  filter(Family != "Empididae") %>%
  filter(Family != "Muscidae") %>%
  filter(Family != "Chironomidae") %>%
  filter(Family != "Dolichopodidae") %>%
  filter(Family != "Scathophagidae") %>%
  filter(Family != "Ceratopogonidae") %>%
  filter(Family != "Hybotidae") %>%
  filter(Family != "Sciaridae") %>%
  filter(Family != "Sphaeroceridae") %>%
  filter(Family != "unknown") %>%
  select(Genus) %>%
  count(Genus, name = "countgenera")

ggplot(KBIMP_Domcontigs_bitinggenus, aes(y=countgenera, x=Genus)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 
#2 mosquito genera and 4 black fly genera

rm(KBIMP_Domcontigs_bitinggenus)

#looking at thr black fly species 

KBIMP_Domcontigs_blackflyspecies <- KBIMP_Domcontigs %>%
  filter(Family == "Simuliidae") %>%
  select(Species) %>%
  count(Species, name = "countspecies")

length(KBIMP_Domcontigs_blackflyspecies$Species) #12 species of black fly are seen here but 2 of them are unknowns 

ggplot(KBIMP_Domcontigs_blackflyspecies, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

rm(KBIMP_Domcontigs_blackflyspecies)

#looking at the mosquito species 

KBIMP_Domcontigs_mosquitoes <- KBIMP_Domcontigs %>%
  filter(Family == "Culicidae") %>%
  select(Species) %>%
  count(Species, name = "countspecies")

length(KBIMP_Domcontigs_mosquitoes$Species) #2 species of black fly are seen here but most of them are unknowns 

ggplot(KBIMP_Domcontigs_mosquitoes, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

rm(KBIMP_Domcontigs_mosquitoes)


#### - Mosquitoes - Phylogenetic analysis and resolution of undefined species ---- 

#filtering for Culicidae, and defining the root as a chirmoidae sequence in the data set 

kbimp_mos_DNA_df <- KBIMP_Domcontigs %>%
  mutate(Family = if_else(Sample == "KBIMP-_005_A2", "Root", Family)) %>%
  filter(Family %in% c("Culicidae", "Root")) %>%
  mutate(Sample = if_else(Sample == "KBIMP-_005_A2", "Root", Sample)) %>%
  select(Sample, Sequence) 


view(kbimp_mos_DNA_df)

# Convert the 'sequence' column to a DNAStringSet
kbimp_mos_DNA <- DNAStringSet(kbimp_mos_DNA_df$Sequence)

# Optionally, assign names to the DNAStringSet
names(kbimp_mos_DNA) <- kbimp_mos_DNA_df$Sample

#preforming an alignment 

alighned_kbimpmos_DNA <- DNAStringSet(muscle::muscle(kbimp_mos_DNA))

#converting fasta file into a phyDat format 

kbimp_mos_phydat <- as.phyDat(alighned_kbimpmos_DNA, type = "DNA")
class(kbimp_mos_phydat) # is a "phyDat" object
length(kbimp_mos_phydat) # has the  species as seen before

#determining the best model for the ML tree 

#modelTest <- modelTest(KGLGTK_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) # JC is the best model because it has the lowest AIC, second lowest BIC and high loglik


# making a tree with maximum likelihood using the JC model
# creating a distance matrix
dist.kbimp.mos <- dist.ml(kbimp_mos_phydat, ratio = TRUE, model = "JC69") 

# Hierarchical clustering
hc <- hclust(dist.kbimp.mos)

# function to find medoid sequence in a cluster
find_medoid <- function(distmat, members) {
  mat <- as.matrix(distmat)
  members <- intersect(members, rownames(mat)) # only keep matching names
  submat <- mat[members, members, drop = FALSE]
  sums <- rowSums(submat)
  medoid <- names(which.min(sums))
  return(medoid)
}


# get cluster memberships
groups <- cutree(hc, k = 14) # I did this mulitple different times to determine the optimal number of groups based on my data 
table(groups) #this displays the number of sequences in each group
cluster_sizes <- table(groups) #this puts that information into a 
# find medoids for each cluster
medoids <- sapply(unique(groups), function(g) {
  members <- names(groups[groups == g])
  find_medoid(dist.kbimp.mos, members)
})

medoids


# convert subset to a phyDat object
medoid_phy <- kbimp_mos_phydat[medoids]

# check that this worked 
medoid_phy #14 sequences with 659 character and 152 different site patterns.

#create a new dist matrix
dist.medoid.mos <- dist.ml(medoid_phy, ratio = TRUE, model = "JC69") 

# creating a tree using the neighbor joining method
NJtree.kbimp.mos <- NJ(dist.medoid.mos)
plot(NJtree.kbimp.mos) #a tree was produced
length(NJtree.kbimp.mos$tip.label) #1272 tips as expected 

# Fit the initial tree using a simple pml 
pml.tree.kbimp.mos <- pml(NJtree.kbimp.mos, kbimp_mos_phydat, k = 4, model = "GTR+I", method = "unrooted")
plot(pml.tree.kbimp.mos$tree) #a tree was produced
pml.tree.kbimp.mosop <- pml_bb(pml.tree.kbimp.mos, model = "GTR+I")

#bootstrapping to determine relability and to take the best tree
bs <- bootstrap.pml(pml.tree.kbimp.mosop, bs = 1000, optNni = TRUE, multicore = TRUE)
tree_with_bs <- plotBS(pml.tree.kbimp.mos$tree, bs)

#rooting the tree
rooted.bstree.mos <- root(tree_with_bs, outgroup = "Root", resolve.root = TRUE)
plot(rooted.bstree.mos)
tree_with_bs <- plotBS(rooted.bstree.mos, bs)

#After BLASTing these sequences I change the name of these 14 sequences and add the species name 

speciesdata <- KBIMP_Domcontigs %>%
  filter(Sample %in% c("KBIMP_004_F7" ,  "KBIMP_002_F10" , "KBIMP_004_B1"  , "KBIMP_003_G5" ,
                       "KBIMP_004_B8"  ,   "KBIMP_002_H11" ,
                        "KBIMP-_005_A2"     ,      "KBIMP-_005_D8" ,
                        "KBIMP-_007_G10", "KBIMP-_007_F10" ,"KBIMP-_007_F2", 
                        "KBIMP-_007_H4" , "KBIMP-_007_H11", "KBIMP-_008_A1")) %>%
  select(Sample, Species) %>%
  mutate(Species = if_else(Sample == "KBIMP-_005_A2", "Root", Species)) %>%
  mutate(Sample = if_else(Sample == "KBIMP-_005_A2", "Root", Sample)) %>%
  mutate(Species = if_else(Sample == "KBIMP_002_F10", "Aedes nigripes", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_004_F7", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_004_B8", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_002_H11", "Aedes impiger", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_005_D8", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_C6", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_004_B1", "Aedes nigripes", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_002_G7", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_G10", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_F2", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_A1", "Aedes excrucians", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_H4", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_H11", "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_003_G5", "Aedes punctor/Aedes hexodontus", Species))
  

#now I produce the finalized tree with the species for each group shown on the tips

# the next few steps are done to make sure that the trait data and the tree have the same order of species for the tree construction
tree_tips_mos <- rooted.bstree.mos$tip.label 
sample_name_mos <- speciesdata$Sample # getting the sample names from the trait data
# checking if all species names are present in the tree
all(tree_tips_mos %in% sample_name_mos) #FALSE not all sample_names in tree tips
all(sample_name_mos %in% tree_tips_mos) #TRUE all the sample names are in the three tips
setdiff(tree_tips_mos, sample_name_mos)
setdiff(sample_name_mos, tree_tips_mos)

speciesdata <- speciesdata[match(tree_tips_mos, speciesdata$Sample), ] %>%
  rownames_to_column(var = "cluster")

speciesdata$cluster_size <- cluster_sizes[as.character(speciesdata$cluster)]

# Assign colors to groups
group_colors <- c("Aedes nigripes"="#00D400","Aedes punctor/Aedes hexodontus"= "#FF2A7F","Aedes impiger"= "#008066","Aedes excrucians"="#FF6600", "Culiseta alaskaensis"="#0069FF","Root"= "#C00000" )
# Add colors to the trait data
speciesdata$color <- group_colors[speciesdata$Species]

#plotting the tree with trait data 
node <- 1:Ntip(rooted.bstree.mos)


#produce the tree
mosquitotree2024 <- (
  ggtree(rooted.bstree.mos, layout = "dendrogram", branch.length = TRUE) +
    geom_tippoint(aes(color = speciesdata$Species[node]), size = 2) +
    scale_color_manual(values = group_colors, name = "Species") +
    geom_text(aes(label = speciesdata$cluster_size[node]), vjust = 1.7, size =3) +
    theme_tree2(legend.position = "right") +
    theme(
      legend.title = element_text(size = 12),      
      legend.text = element_text(size = 10)
    )
)

mosquitotree2024

#save the tree
ggsave("../plots/mosphylomoregroups.png", plot = mosquitotree2024, width = 6, height = 6, dpi = 300)

#### ---- Mosquitoes - Getting the sequences for haplotype networks ---- 


# Hierarchical clustering
hc <- hclust(dist.kbimp.mos)
groups_haplotype <- cutree(hc, k = 5) #this clustering allows me to take the various subspecies/haplotypes shown on the tree together to be saved as the required file 

# View assignments
table(groups_haplotype)

seqs_punctor_mos <- alighned_kbimpmos_DNA[names(groups_haplotype)[groups_haplotype == 1]]
seqs_nigripes_mos <- alighned_kbimpmos_DNA[names(groups_haplotype)[groups_haplotype == 2]]
seqs_excrucians_mos <- alighned_kbimpmos_DNA[names(groups_haplotype)[groups_haplotype == 4]]


kbimp_mosmetadata <- kbimp2024_sampledata %>%
  inner_join(KBIMP_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, SamplingProtocol, FieldID, ExactSite, Lat, Lon, Family, Sector) %>%
  filter(Family == "Culicidae") %>%
  add_row(SampleID = "Root", Sector = "Root")

kbimp_mosmetadata_punctor <- kbimp_mosmetadata %>%
  filter(SampleID %in% names(seqs_punctor_mos)) %>%
  select(SampleID, ExactSite) %>%
  column_to_rownames(var = "SampleID")

kbimp_mosmetadata_nigirpes <- kbimp_mosmetadata %>%
  filter(SampleID %in% names(seqs_nigripes_mos)) %>%
  select(SampleID, ExactSite) %>%
  column_to_rownames(var = "SampleID")

kbimp_mosmetadata_excrucians <- kbimp_mosmetadata %>%
  filter(SampleID %in% names(seqs_excrucians_mos)) %>%
  select(SampleID, ExactSite) %>%
  column_to_rownames(var = "SampleID")

# Convert to binary matrix
binary_matrix_trait_punctor <- model.matrix(~ ExactSite - 1, data = kbimp_mosmetadata_punctor)
binary_matrix_trait_nigirpes <- model.matrix(~ ExactSite - 1, data = kbimp_mosmetadata_nigirpes)
binary_matrix_trait_excrucians <- model.matrix(~ ExactSite - 1, data = kbimp_mosmetadata_excrucians)

binary_matrix_trait_punctor <- binary_matrix_trait_punctor %>% as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  select("SampleID","ExactSitedew line rd" , "ExactSitefreshwater river","ExactSiteGrenier lake",
         "ExactSitelong pt creek", "ExactSitegravel pit", "ExactSite4 mile bay",
         "ExactSitebehind town", "ExactSiteother side of heart lake", "ExactSitetrail to 4 mile")

binary_matrix_trait_nigirpes <- binary_matrix_trait_nigirpes %>% as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  select("SampleID","ExactSitedew line rd" , "ExactSitefreshwater river","ExactSiteGrenier lake",
         "ExactSitelong pt creek",  "ExactSitetrail to 4 mile")

binary_matrix_trait_excrucians <- binary_matrix_trait_excrucians %>% as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  select("SampleID","ExactSite4 mile bay" , "ExactSitebehind town","ExactSiteother side of heart lake",
         "ExactSitetrail to 4 mile")


#writing the binary matrices and nexus files 
write_tsv(binary_matrix_trait_punctor, "../processed-data/mosmetadata_punctor.tsv")
write.nexus.data(seqs_punctor_mos, file = "../processed-data/output_alignment_punctor.nex")
write_tsv(binary_matrix_trait_nigirpes, "../processed-data/mosmetadatag_nigirpes.tsv")
write.nexus.data(seqs_nigripes_mos, file = "../processed-data/output_alignment_nigirpes.nex")
write_tsv(binary_matrix_trait_excrucians, "../processed-data/mosmetadataexcrucians.tsv")
write.nexus.data(seqs_excrucians_mos, file = "../processed-data/output_alignment_excrucians.nex")


#### - Black flies - Phylogenetic analysis and resolution of undefined species ----

#filtering for Simuliidae, and defining the root as a chirmoidae sequence in the data set 

kbimp_bf_DNA_df <- KBIMP_Domcontigs %>%
  mutate(Family = if_else(Sample == "KBIMP-_005_A2", "Root", Family)) %>%
  filter(Family %in% c("Simuliidae", "Root")) %>%
  mutate(Sample = if_else(Sample == "KBIMP-_005_A2", "Root", Sample)) %>%
  select(Sample, Sequence) 
 

view(kbimp_bf_DNA_df)

# Convert the 'sequence' column to a DNAStringSet
kbimp_bf_DNA <- DNAStringSet(kbimp_bf_DNA_df$Sequence)

# Optionally, assign names to the DNAStringSet
names(kbimp_bf_DNA) <- kbimp_bf_DNA_df$Sample

#preforming an alignment 

alighned_kbimpbf_DNA <- DNAStringSet(muscle::muscle(kbimp_bf_DNA))

#converting fasta file into a phyDat format 

kbimp_bf_phydat <- as.phyDat(alighned_kbimpbf_DNA, type = "DNA")
class(kbimp_bf_phydat) # is a "phyDat" object
length(kbimp_bf_phydat) # has the  species as seen before

#determining the best model for the ML tree 

#modelTest <- modelTest(KGLGTK_mos_phydat, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR")) # JC is the best model because it has the lowest AIC, second lowest BIC and high loglik


# making a tree with maximum likelihood using the JC model
# creating a distance matrix
dist.kbimp.bf <- dist.ml(kbimp_bf_phydat, ratio = TRUE, model = "JC69") 

# Hierarchical clustering
hc_bf <- hclust(dist.kbimp.bf)

# get cluster memberships
groups <- cutree(hc_bf, k = 20) # I did this multiple different times to determine the optimal number of groups based on my data 
table(groups) #this displays the number of sequences in each group
cluster_sizes <- table(groups) #this puts that information into a 
# find medoids for each cluster
medoids <- sapply(unique(groups), function(g) {
  members <- names(groups[groups == g])
  find_medoid(dist.kbimp.bf, members)
})

medoids


# convert subset to a phyDat object
medoid_phy <- kbimp_bf_phydat[medoids]

# check that this worked 
medoid_phy #20 sequences with 662 character and 252 different site patterns.

#create a new dist matrix
dist.medoid.bf <- dist.ml(medoid_phy, ratio = TRUE, model = "JC69") 

# creating a tree using the neighbor joining method
NJtree.kbimp.bf <- NJ(dist.medoid.bf)
plot(NJtree.kbimp.bf) #a tree was produced
length(NJtree.kbimp.bf$tip.label) #1272 tips as expected 

# Fit the initial tree using a simple pml 
pml.tree.kbimp.bf <- pml(NJtree.kbimp.bf, kbimp_bf_phydat, k = 4, model = "GTR+I", method = "unrooted")
plot(pml.tree.kbimp.bf$tree) #a tree was produced
pml.tree.kbimp.bfop <- pml_bb(pml.tree.kbimp.bf, model = "GTR+I")

bs <- bootstrap.pml(pml.tree.kbimp.bfop, bs = 1000, optNni = TRUE, multicore = TRUE)
tree_with_bsbf <- plotBS(pml.tree.kbimp.bf$tree, bs)

#rooting the tree
rooted.bstree.bf <- root(tree_with_bsbf, outgroup = "Root", resolve.root = TRUE)
plot(rooted.bstree.bf)
tree_with_bs <- plotBS(rooted.bstree.bf, bs)

#After BLASTing these sequences I change the name of these 14 sequences and add the species name 

speciesdata <- KBIMP_Domcontigs %>%
  filter(Sample %in% c("KBIMP_01_B8" ,  "KBIMP_002_D12" , "KBIMP_002_E5"  , "KBIMP-_005_F4" ,
                       "KBIMP-_006_D8"  ,   "KBIMP-_007_E10" ,
                       "KBIMP-_006_C5"     ,      "KBIMP_003_B12" ,
                       "KBIMP_003_E7", "KBIMP_003_F2" ,"KBIMP-_005_A2",
                       "KBIMP-_007_C11" , "KBIMP-_007_C5", "KBIMP-_007_D4", "KBIMP-_007_E9",
                       "KBIMP-_008_A3", "KBIMP-_008_B2", "KBIMP-_008_D10", "KBIMP-_008_D4", "KBIMP-_008_D5")) %>%
  select(Sample, Species) %>%
  mutate(Species = if_else(Sample == "KBIMP-_005_A2", "Root", Species)) %>%
  mutate(Sample = if_else(Sample == "KBIMP-_005_A2", "Root", Sample)) %>%
  mutate(Species = if_else(Sample == "KBIMP_002_E5", "Simulium undescribed", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_006_D8", "Simulium noelleri", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_E10", "Simulium tuberosum complex", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP_003_F2", "Simulium decorum", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_C11", "Simulium arcticum complex", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_007_E9", "Simulium tuberosum complex", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_A3", "Simulium venustum complex", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_B2", "Simulium subpusillum", Species)) %>%
  mutate(Species = if_else(Sample == "KBIMP-_008_D4", "Simulium bicorne", Species)) 

#making a tree with the sampling locations labels on the tree tips 

tree_tips_bf <- rooted.bstree.bf$tip.label 
sample_name_bf <- speciesdata$Sample # getting the sample names from the trait data
# checking if all species names are present in the tree
# checking if all species names are present in the tree
all(tree_tips_bf %in% sample_name_bf) #FALSE not all sample_names in tree tips
all(sample_name_bf %in% tree_tips_bf) #TRUE all the sample names are in the three tips
setdiff(tree_tips_bf, sample_name_bf)
setdiff(sample_name_bf, tree_tips_bf)

speciesdata <- speciesdata[match(tree_tips_bf, speciesdata$Sample), ] %>%
  rownames_to_column(var = "cluster")

speciesdata$cluster_size <- cluster_sizes[as.character(speciesdata$cluster)]

# Assign colors to groups
group_colors <- c("Metacnephia borealis"= "chartreuse2", 
                  "Cnephia eremites" = "#FF2A7F", 
                  "Root" = "#C00000",
                  "Simulium noelleri" = "cadetblue4", 
                  "Simulium decimatum" = "gold", 
                  "Simulium baffinense" = "orangered", 
                  "Simulium subpusillum" = "blueviolet", 
                  "Simulium congareenarum" = "forestgreen",
                  "Simulium malyschevi" = "blue3", 
                  "Simulium undescribed" = "maroon", 
                  "Simulium craigi" = "darkorange2", 
                  "Simulium tuberosum complex"="dodgerblue", 
                  "Stegopterna emergens" = "olivedrab3", 
                  "Simulium bicorne" = "orange4",
                  "Simulium venustum complex"="purple4", 
                  "Simulium decorum" ="coral",
                  "Simulium arcticum complex"= "lightyellow4")
# Add colors to the trait data
speciesdata$color <- group_colors[speciesdata$Species]

#plotting the tree with trait data 
node <- 1:Ntip(rooted.bstree.bf)




bftree2024 <- (
  ggtree(rooted.bstree.bf, layout = "dendrogram", branch.length = TRUE) +
    geom_tiplab(size =2) +
    geom_tippoint(aes(color = speciesdata$Species[node]), size = 2) +
    scale_color_manual(values = group_colors, name = "Species") +
    geom_text(aes(label = speciesdata$cluster_size[node]), vjust = 1.7, size =3) +
    theme_tree2(legend.position = "right") +
    theme(
      legend.title = element_text(size = 12),      
      legend.text = element_text(size = 10)
    )
)
bftree2024
ggsave("../plots/bfphylo.png", plot = bftree2024, width = 6, height = 6, dpi = 300)


#### - Black flies - Getting the sequences for haplotype networks ---- 

# Hierarchical clustering
hc <- hclust(dist.kbimp.bf)
groups <- cutree(hc, k = 20)

# View assignments
table(groups)

seqs_group1 <- alighned_kbimpbf_DNA[names(groups)[groups == 1]] #174: Metacnephia borealis - same with 16
seqs_group2 <- alighned_kbimpbf_DNA[names(groups)[groups == 2]]  #21: Cnephia eremites - same with 16

seqs_group4 <- alighned_kbimpbf_DNA[names(groups)[groups == 4]] #6: Simulium subpusillum - same with 16
seqs_group5 <- alighned_kbimpbf_DNA[names(groups)[groups == 5]] #13: Simulium noelleri - 12 with 16 
seqs_tuberosum1 <- alighned_kbimpbf_DNA[names(groups)[groups == 7]] 
seqs_group8 <- alighned_kbimpbf_DNA[names(groups)[groups == 8]] #1: Stegopterna emergens - same with 16
seqs_group9 <- alighned_kbimpbf_DNA[names(groups)[groups == 9]] #13: Simulium decimatum - same with 16
seqs_group10 <- alighned_kbimpbf_DNA[names(groups)[groups == 10]] #1: Simulium baffinense - same with 16
seqs_decorum  <- alighned_kbimpbf_DNA[names(groups)[groups == 11]] #1: Simulium decorum 
seqs_arcticum1 <- alighned_kbimpbf_DNA[names(groups)[groups == 12]] #6:  Simulium arcticum complex

seqs_group13 <- alighned_kbimpbf_DNA[names(groups)[groups == 13]] #18:Simulium congareenarum
seqs_group14 <- alighned_kbimpbf_DNA[names(groups)[groups == 14]] #11: Simulium malyschevi
seqs_tuberosum2 <- alighned_kbimpbf_DNA[names(groups)[groups == 15]] #3: S tuberosum complex - same with 16
seqs_venustum <- alighned_kbimpbf_DNA[names(groups)[groups == 16]] #2: Simulium venustum complex
seqs_subpusillum <- alighned_kbimpbf_DNA[names(groups)[groups == 17]] #2: Simulium subpusillum
seqs_bicorne <- alighned_kbimpbf_DNA[names(groups)[groups == 19]] #2: Simulium bicorne


seqs_Simulium_undescribed  <- alighned_kbimpbf_DNA[names(groups)[groups == 3]] #86: Simulium sp. Unclassified - same with 16
seqs_noelleri <- alighned_kbimpbf_DNA[names(groups)[groups == 5]] #13: Simulium noelleri - 12 with 16 
#with 15 cuts it puts the norelli and decorum together with an addtional cut it sperates them
#with 17 cuts it seperates the tubersum complex into a group of 8 and a group of three - the group of 6 is" Simulium annulitarse and the group of three is:Simulium
#tuberosum complex and unspecified 


#### - Ecological analysis ----



KBIMP_updatedspecies <- KBIMP_Domcontigs %>%
  mutate(Species = ifelse(Sample %in%names(seqs_punctor_mos) , "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_nigripes_mos) , "Aedes impiger/Aedes nigripes", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_excrucians_mos) , "Aedes excrucians", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_Simulium_undescribed) , "Simulium undescribed", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_noelleri) , "Simulium noelleri", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_tuberosum1) , "Simulium tuberosum complex", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_tuberosum2) , "Simulium tuberosum complex", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_decorum) , "Simulium decorum", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_arcticum1) , "Simulium arcticum complex", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_venustum) , "Simulium venustum complex", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_subpusillum) , "Simulium subpusillum", Species)) %>%
  mutate(Species = ifelse(Sample %in%names(seqs_bicorne) , "Simulium bicorne", Species)) 

#looking at thr black fly species 

KBIMP_Domcontigs_blackflyspecies <- KBIMP_updatedspecies %>%
  filter(Family == "Simuliidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP_Domcontigs_blackflyspecies$Species) #12 species of black fly are seen here but 2 of them are unknowns 

ggplot(KBIMP_Domcontigs_blackflyspecies, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

rm(KBIMP_Domcontigs_blackflyspecies)

#looking at the mosquito species 

KBIMP_Domcontigs_mosquitoes <- KBIMP_updatedspecies %>%
  filter(Family == "Culicidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP_Domcontigs_mosquitoes$Species) #2 species of black fly are seen here but most of them are unknowns 

ggplot(KBIMP_Domcontigs_mosquitoes, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

rm(KBIMP_Domcontigs_mosquitoes)


#now I make this data into a species by site matrix with each repitition at each site combined

species_matrix_CBAY_sites <- KBIMP_updatedspecies %>%
  full_join(kbimp2024_sampledata, join_by(Sample == SampleID)) %>%
  filter(Family %in% c("Culicidae", "Simuliidae")) %>%
  select(ExactSite, Species, Sector) %>%
  filter(Sector == "CBAY") %>%
  dplyr::count(ExactSite, Species, name = "Abundance") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "ExactSite") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))


#run beatdiverity analysis 

betacbaybf <- beta.pair(species_matrix_CBAY_sites, index.family = "sorensen")  


#converting this P/A into the long format

mat_long <- species_matrix_CBAY_sites %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site") %>%
  pivot_longer(-Site, names_to = "Species", values_to = "Presence") %>%
  mutate(Presence = factor(Presence, levels = c(0,1)))

#ploting the tile plot

tileplot <-ggplot(mat_long, aes(x = Site, y = Species, fill = Presence)) +
  scale_x_discrete(labels = c("1st Lake"="FIL", "dew line rd"="DEW", 
                              "freshwater river"="FWR", "Grenier lake"="GNL", 
                              "long pt creek"="LPC", "gravel pit" = "GVP")) +
  geom_tile(color = "black", size = 0.75) +
  scale_fill_manual(values = c("0" = "#0072B2", "1" = "#D55E00")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black", margin = margin(t = 5)), axis.title.x = element_text(size = 12),      # X-axis title font size
        axis.title.y = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        legend.title = element_text(size = 12),      
        legend.text = element_text(size = 10)) +
  labs(x = "Sampling Location", y = "Species")

tileplot

ggsave("../plots/tileplotsl.png", plot = tileplot , width = 6, height = 3, dpi = 300)

#now I do this again but for KGLKTK 

species_matrix_KGLTK_sites <- KBIMP_updatedspecies %>%
  full_join(kbimp2024_sampledata, join_by(Sample == SampleID)) %>%
  filter(Family %in% c("Culicidae", "Simuliidae")) %>%
  select(ExactSite, Species, Sector) %>%
  filter(Sector == "KGLTK") %>%
  dplyr::count(ExactSite, Species, name = "Abundance") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "ExactSite") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))


#run beatdiverity analysis 

betacbayKGLTK <- beta.pair(species_matrix_KGLTK_sites, index.family = "sorensen")  


#converting this P/A into the long format

mat_long <- species_matrix_KGLTK_sites %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site") %>%
  pivot_longer(-Site, names_to = "Species", values_to = "Presence") %>%
  mutate(Presence = factor(Presence, levels = c(0,1)))

#ploting the tile plot

tileplot <-ggplot(mat_long, aes(x = Site, y = Species, fill = Presence)) +
  scale_x_discrete(labels = c("4 mile bay"="4MB", "behind town"="BHT", 
                              "other side of heart lake"="OHL", "trail to 4 mile"="T4M", 
                              "behind horseshoe"="BHS")) +
  geom_tile(color = "black", size = 0.75) +
  scale_fill_manual(values = c("0" = "#0072B2", "1" = "#D55E00")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black", margin = margin(t = 5)), axis.title.x = element_text(size = 12),      # X-axis title font size
        axis.title.y = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        legend.title = element_text(size = 12),      
        legend.text = element_text(size = 10)) +
  labs(x = "Sampling Location", y = "Species")

tileplot

ggsave("../plots/tileplotsl.png", plot = tileplot , width = 6, height = 3, dpi = 300)


#now I do betadiverity comparing the two communities

species_matrix_communties <- KBIMP_updatedspecies %>%
  full_join(kbimp2024_sampledata, join_by(Sample == SampleID)) %>%
  filter(Family %in% c("Culicidae", "Simuliidae")) %>%
  select(Species, Sector) %>%
  dplyr::count(Sector, Species, name = "Abundance") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sector") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))


#run beatdiverity analysis 

betacommunties <- beta.pair(species_matrix_communties, index.family = "sorensen")  


#converting this P/A into the long format

mat_long <- species_matrix_KGLTK_sites %>%
  as.data.frame() %>%
  rownames_to_column(var = "Site") %>%
  pivot_longer(-Site, names_to = "Species", values_to = "Presence") %>%
  mutate(Presence = factor(Presence, levels = c(0,1)))

#ploting the tile plot

tileplot <-ggplot(mat_long, aes(x = Site, y = Species, fill = Presence)) +
  scale_x_discrete(labels = c("4 mile bay"="4MB", "behind town"="BHT", 
                              "other side of heart lake"="OHL", "trail to 4 mile"="T4M", 
                              "behind horseshoe"="BHS")) +
  geom_tile(color = "black", size = 0.75) +
  scale_fill_manual(values = c("0" = "#0072B2", "1" = "#D55E00")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black", margin = margin(t = 5)), axis.title.x = element_text(size = 12),      # X-axis title font size
        axis.title.y = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        legend.title = element_text(size = 12),      
        legend.text = element_text(size = 10)) +
  labs(x = "Sampling Location", y = "Species")

tileplot

ggsave("../plots/tileplotsl.png", plot = tileplot , width = 6, height = 3, dpi = 300)





