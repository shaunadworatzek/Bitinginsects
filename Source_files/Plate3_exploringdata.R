#### - Cleaning the data ----


#load required packages

library(readr)
library(tidyverse)

#set working directory to source file location

setwd("C:/Users/sdwor/OneDrive - University of Guelph/1- MASTERS DEGREE/Bitinginsects/Source_files")

#open data into R

Plate3_Allcontigs <- read_tsv(file = "../Data/CBAYPlate3_seq3_100folddil/Shauna_100dilution/Shauna_100dilution_TaxonomicAssignments_AllContigs.tsv")

View(Plate3_Allcontigs)

#converting the data into the long format 

Plate3_Allcontigs_2 <- Plate3_Allcontigs %>%
  mutate(Genus_Species = str_c(Genus, Species, sep = "_")) %>%
  select(Sample, Genus_Species, ReadCount) %>%
  group_by(Sample, Genus_Species) %>%
  summarise(
    Genus_Species = str_c(unique(Genus_Species), collapse = ", "),  # Combine species names
    Total_Reads = sum(ReadCount, na.rm = TRUE)  # Sum the read counts
  ) %>%
  ungroup()  %>%
  group_by(Sample) %>%
  filter(n() > 1) %>%
  ungroup()  #Some of the samples which had more than one contig were checked via the pictures online -- almost all seemed to correspond to the one with the highest number of reads 
view(Plate3_Allcontigs_2)

#now looking at the diveristy from broader groups to see if we got any interactions:

#phylum
Plate3_Allcontigs_phylum <- Plate3_Allcontigs %>%
  select(Phylum, Class, Order) %>%
  count(Phylum, name = "countphylum")
view(Plate3_Allcontigs_phylum)

ggplot(Plate3_Allcontigs_phylum, aes(y=countphylum, x=Phylum)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 

#class
Plate3_Allcontigs_class <- Plate3_Allcontigs %>%
  select(Phylum, Class, Order) %>%
  count(Class, name = "countclass")
view(Plate3_Allcontigs_class)

ggplot(Plate3_Allcontigs_class, aes(y=countclass, x=Class)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 

#order
Plate3_Allcontigs_order <- Plate3_Allcontigs %>%
  select(Phylum, Class, Order) %>%
  count(Order, name = "countorder")
view(Plate3_Allcontigs_order)

ggplot(Plate3_Allcontigs_order, aes(y=countorder, x=Order)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 

#The only possible coinhabitiing species were the bacteria:Rickettsiales (although this is exciting because this group includes wolbacia this was on a muscidae so not super interesting)
#we also had a Trombidiformes (type of mite) specificaaly the species it identified was a water mite


#now lets open the primary contigs to see what species we actually had in each well 

Plate3_Domcontigs <- read_tsv(file = "../Data/CBAYPlate3_seq3_100folddil/Shauna_100dilution/Shauna_100dilution_TaxonomicAssignments_DominantContigs.tsv")
view(Plate3_Domcontigs)


#lets look at the amount of wrong IDs

#looking at it by Family
Plate3_Domcontigs_family <- Plate3_Domcontigs %>%
  select(Family) %>%
  count(Family, name = "countfamily")
view(Plate3_Domcontigs_family)

ggplot(Plate3_Domcontigs_family, aes(y=countfamily, x=Family)) +
  geom_col(fill = "skyblue3") +
  theme_bw() #we did a good job only 7 which werent mosquitoes or black flies 


#now lets look at the genus diveristy of only black flies and mosquitoes

Plate3_Domcontigs_bitinggenus <- Plate3_Domcontigs %>%
  filter(Family != "Dolichopodidae") %>%
  filter(Family != "Muscidae") %>%
  filter(Family != "Scathophagidae") %>%
  select(Genus) %>%
  count(Genus, name = "countgenera")

ggplot(Plate3_Domcontigs_bitinggenus, aes(y=countgenera, x=Genus)) +
  geom_col(fill = "skyblue3") +
  theme_bw() #we did a good job only 7 which werent mosquitoes or black flies 
#all mosqutioes were aedes - but we saw a split in black flies between Metacnephia and Simulium 

#it seems like all the mosquitoes were only identified to the genus level so lets look at the species of black flies

Plate3_Domcontigs_blackflyspecies <- Plate3_Domcontigs %>%
  filter(Family == "Simuliidae") %>%
  select(Species) %>%
  count(Species, name = "countspecies")

ggplot(Plate3_Domcontigs_blackflyspecies, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw()




