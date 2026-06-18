
#Data required to run analysis is the megtabarcoding data file from 2025:
#KBIMP2025_insectCOI_OTUDetails.tsv, and the two metadata files from PCR 
#and extraction to give information on controls extractiondata2025.csv and 
#PCRcontrol_data.csv. 

#Output will be the filtered metabarcoding file which has the controls filtered 
#and removed

library(stringr)
library(tidyverse)
library(readr)
library(viridis)
library(ggplot2)   

#### - PART 1 - 2024 barcoding data  ----

##### Filtering negative controls #####

#input raw sequence data 

KGLKTK2024_COI <- read_tsv(file = "raw-data2/KGLKTK_2024_OTUDetails.tsv")
CBAY2024_plates12456 <- read_tsv(file = "raw-data2/CBAY2024_AllPlates_OTUDetails.tsv")
CBAY2024_plate3 <- read_tsv(file = "raw-data2/Shauna_CBAY2024_Plate3_OTUDetails.tsv")
problemsamples <- read_csv(file = "raw-data2/problemsamples.csv")

length(KGLKTK2024_COI$Sample) #232
length(CBAY2024_plates12456$Sample) #635
length(CBAY2024_plate3$Sample) #261

#Filtering out/ accounting for negative controls 
#negative control on plate 3

CBAY2024_plate3_fil <- CBAY2024_plate3 %>%
  group_by(Species, Genus) %>%
  mutate(
    control_reads = sum(Read_Count[Sample == "CONTROL_05_H12"], na.rm = TRUE),
    Read_Count = Read_Count - control_reads
  ) %>%
  filter(!Read_Count <= 3) %>%
  filter(Barcode_Status == "target") %>%
  select(Sample, Read_Count, Kingdom, Phylum, Class, Order, Family, Genus, 
         Species, Probability_Genus, Probability_Species, Sequence)
#filter(!Sample %in% c("KBIMP_003_G3")) #must have been a reason I filtered this out have to figure that out???

#negative controls on plate 1,2,4,5,6

CBAY2024_plates12456_fil <- CBAY2024_plates12456 %>%
  filter(!Sample %in% c("KBIMP_004_H11", "KBIMP_004_D12", "KBIMP_004_E12", "KBIMP-_006_G4",
                        "KBIMP_004_F12", "KBIMP_004_G12", "KBIMP_004_H10", "KBIMP-_006_H1",
                        "KBIMP-_006_H10", "KBIMP-_006_H4", "KBIMP-_006_H5", "KBIMP-_006_H6",
                        "KBIMP-_006_H8", "KBIMP_01_G11", "KBIMP_01_G12", "KBIMP_01_H1", "KBIMP_01_G11", "KBIMP-_006_H11", "KBIMP_002_B11")) %>%
  filter(!Sample %in% problemsamples$`problem samples`) %>%
  mutate(`Plate #` = case_when(grepl("^KBIMP", Sample) ~ sub("_[A-H][0-9]{1,2}$", "", Sample),
                               
                               grepl("^CONTROL_", Sample) ~paste0("KBIMP-_",str_pad(as.integer(str_extract(Sample, "[0-9]+$")),3,
                                                                                    pad = "0")), TRUE ~ NA_character_)) %>%
  mutate(`Plate #` = gsub("[-_]", "", `Plate #`)) %>%
  group_by(`Plate #`, Species, Genus) %>%
  mutate(
    control_reads = sum(Read_Count[grepl("^CONTROL", Sample)], na.rm = TRUE),
    Read_Count = Read_Count - control_reads
  ) %>%
  ungroup() %>%
  filter(Read_Count > 3) %>%
  filter(Barcode_Status == "target") %>%
  select(Sample, Read_Count, Kingdom, Phylum, Class, Order, Family, Genus, 
         Species, Probability_Genus, Probability_Species, Sequence)

#negative controls for Kugluktuk plate 

KGLKTK2024_COI_fil <- KGLKTK2024_COI %>%
  filter(!Sample %in% c("KBIMP-_008_F4", "KBIMP-_008_F5", "KBIMP_01_B3", "KBIMP_01_F8", "KBIMP_01_B11",
                        "KBIMP_01_D2", "KBIMP_01_E5", "KBIMP_01_F7", "KBIMP_01_B12")) %>%
  mutate(`Plate #` = case_when(grepl("^KBIMP", Sample) ~ sub("_[A-H][0-9]{1,2}$", "", Sample),
                               Sample == "CONTROL_H12_Plate1" ~ "KBIMP-_007",
                               Sample == "CONTROL_H12_Plate2" ~ "KBIMP-_008", TRUE ~ NA_character_)) %>%
  group_by(`Plate #`, Species, Genus) %>%
  mutate(
    control_reads = sum(Read_Count[grepl("^CONTROL", Sample)], na.rm = TRUE),
    Read_Count = Read_Count - control_reads
  ) %>%
  ungroup() %>%
  filter(Read_Count > 3) %>%
  filter(Barcode_Status == "target") %>%
  select(Sample, Read_Count, Kingdom, Phylum, Class, Order, Family, Genus, 
         Species, Probability_Genus, Probability_Species, Sequence)


length(KGLKTK2024_COI_fil$Sample) #152 lost 48 (some of these were due to sampling issues and not controls)
length(CBAY2024_plates12456_fil$Sample) #359 lost 76 (but this includes the 7 repeats)
length(CBAY2024_plate3_fil$Sample) #95 lost 0 (there were doubles of some samples)

##### Combining three sequencing runs #####

CBAY2024 <- rbind(CBAY2024_plates12456_fil, CBAY2024_plate3_fil)

KBIMP2024 <- rbind(CBAY2024, KGLKTK2024_COI_fil)

write_tsv(KBIMP2024, "processed-data/KBIMP2024_filteredCOI.tsv")

##### Investigating species in data #####

#lets look at the amount of wrong IDs (or IDs that aren't black flies/mosquitoes)

#looking at it by Family
KBIMP2024_family <- KBIMP2024 %>%
  select(Family) %>%
  dplyr::count(Family, name = "countfamily")

ggplot(KBIMP2024_family, aes(y=countfamily, x=Family)) +
  geom_col(fill = "skyblue3") +
  theme_bw()  #7 families which weren't mosquitoes or black flies 

#now lets look at the genus diversity of only black flies and mosquitoes

KBIMP2024_bitinggenus <- KBIMP2024 %>%
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
  dplyr::count(Genus, name = "countgenera")

ggplot(KBIMP2024_bitinggenus, aes(y=countgenera, x=Genus)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 
#2 mosquito genera and 4 black fly genera

#looking at thr black fly species 

KBIMP2024_blackflyspecies <- KBIMP2024 %>%
  filter(Family == "Simuliidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP2024_blackflyspecies$Species) #13 species of black fly are seen here but 2 of them are unknowns 

ggplot(KBIMP2024_blackflyspecies, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

#looking at the mosquito species 

KBIMP2024_mosquitoes <- KBIMP2024 %>%
  filter(Family == "Culicidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP2024_mosquitoes$Species) #4 species of black fly are seen here but most of them are unknowns 

ggplot(KBIMP2024_mosquitoes, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

rm(KBIMP2024_mosquitoes, CBAY2024, problemsamples,
   KBIMP2024_blackflyspecies, KBIMP2024_bitinggenus, KGLKTK2024_COI_fil, 
   CBAY2024_plates12456_fil, CBAY2024_plate3_fil, KGLKTK2024_COI, 
   CBAY2024_plates12456, CBAY2024_plate3, KBIMP2024_family)



#### - PART 2 - Filtering of 2025 metabarcoding data -----

KBIMP2025_COI <- read_tsv(file = "raw-data2/KBIMP2025_insectCOI_OTUDetails.tsv")

##### Filtering for negative controls #####

#opening negative control data into R 

df_extractioncontrols <- read_csv(file = "raw-data2/extractiondata2025.csv")
df_PCRcontrols <- read_csv(file = "raw-data2/PCRcontrol_data.csv")

# filtering out PCR controls 

KBIMP2025_COI_fil <- KBIMP2025_COI %>%
  left_join(df_PCRcontrols, join_by(Sample == Sample), relationship = "many-to-many")  %>%
  group_by(`Plate #`, Species, Genus) %>%
  mutate(
    PCRC_reads = sum(Read_Count[grepl("^PCRC", Sample)], na.rm = TRUE),
    Read_Count = Read_Count - PCRC_reads
  ) %>%
  ungroup() %>%
  filter(!Read_Count <= 3)

#Filtering out Extraction controls
  
KBIMP2025_COI_fil2 <- KBIMP2025_COI_fil %>%
  left_join(df_extractioncontrols , join_by(Sample == Sample), relationship = "many-to-many")  %>%
  group_by(`Plate #`, Species, Genus) %>%
  mutate(
    EXC_reads = sum(Read_Count[grepl("^EXC", Sample)], na.rm = TRUE),
    Read_Count = Read_Count - EXC_reads
  ) %>%
  ungroup() %>%
  filter(!Read_Count <= 3)

length(unique(KBIMP2025_COI_fil2$Sample)) 

###### investigating diversity ######

#lets look at the amount of wrong IDs (or IDs that aren't black flies/mosquitoes)

#looking at it by Family
KBIMP2025_COI_family <- KBIMP2025_COI_fil %>%
  select(Family) %>%
  dplyr::count(Family, name = "countfamily")

ggplot(KBIMP2025_COI_family, aes(y=countfamily, x=Family)) +
  geom_col(fill = "skyblue3") +
  theme_bw()  #7 families which weren't mosquitoes or black flies 

#now lets look at the genus diversity of only black flies and mosquitoes

KBIMP2025_bitinggenus <- KBIMP2025_COI_fil2 %>%
  filter(Family != "Empididae") %>%
  filter(Family != "Muscidae") %>%
  filter(Family != "Chironomidae") %>%
  filter(Family != "Dolichopodidae") %>%
  filter(Family != "Scathophagidae") %>%
  filter(Family != "Ceratopogonidae") %>%
  filter(Family != "Hybotidae") %>%
  filter(Family != "Sciaridae") %>%
  filter(Family != "Sphaeroceridae") %>%
  filter(Family != "Megachilidae") %>%
  filter(Family != "unknown") %>%
  select(Genus) %>%
  dplyr::count(Genus, name = "countgenera")

ggplot(KBIMP2025_bitinggenus, aes(y=countgenera, x=Genus)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 
#2 mosquito genera and 4 black fly genera

#looking at thr black fly species 

KBIMP2025_blackflyspecies <- KBIMP2025_COI_fil2 %>%
  filter(Family == "Simuliidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP2025_blackflyspecies$Species) #16 species of black fly are seen here but 2 of them are unknowns
#4 more species than what was seen in 2024

ggplot(KBIMP2025_blackflyspecies, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

#looking at the mosquito species 

KBIMP2025_mosquitoes <- KBIMP2025_COI_fil2  %>%
  filter(Family == "Culicidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP2025_mosquitoes$Species) #4 species of black fly are seen here but most of them are unknowns 

ggplot(KBIMP2025_mosquitoes, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

write_tsv(KBIMP2025_COI_fil2, "processed-data/KBIMP2025_filteredCOI.tsv")

rm(df_extractioncontrols, df_PCRcontrols, KBIMP2025_COI, KBIMP2025_COI_fil,
   KBIMP2025_mosquitoes, KBIMP2025_blackflyspecies,KBIMP2025_bitinggenus, 
   KBIMP2025_COI_family)





