


#### - Preparing the site metadata for analysis ----

kbimp2024_sampledata <- read_csv(file = "../raw-data/KBIMP2024_specimendata.csv")
kbimp2024_sitesnamesfixed <- read_csv(file = "../raw-data/KBIMP_meta_sitenamesfixed.csv")

view(kbimp2024_sampledata)

KBIMP_2024_siteslook <- kbimp2024_sampledata %>%
  group_by(Site) %>%
  summarise(countsites = n_distinct(FieldID)) 


ggplot(KBIMP_2024_siteslook, aes(y= countsites, x= Site)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(KBIMP_2024_siteslook)

#combine good sites sheet with other sampledata sheet

kbimp2024_sampledata <- kbimp2024_sitesnamesfixed %>%
  inner_join(kbimp2024_sampledata, join_by(FieldID == FieldID))


KBIMP_2024_siteslook <- kbimp2024_sampledata %>%
  group_by(ExactSite) %>%
  summarise(countsites = n_distinct(FieldID)) 


ggplot(KBIMP_2024_siteslook, aes(y= countsites, x= ExactSite)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(KBIMP_2024_siteslook)

#looking at the collectionmethods 

KBIMP_2024_methodlook <- kbimp2024_sampledata %>%
  group_by(SamplingMethod) %>%
  summarise(countsites = n_distinct(FieldID)) 


ggplot(KBIMP_2024_methodlook, aes(y= countsites, x= SamplingMethod)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(KBIMP_2024_methodlook)

#fixing the samples which were switched during extractions

kbimp2024_sampledata <- kbimp2024_sampledata %>%
  mutate(SampleID= gsub("KBIMP_004_H11", "KBIMP-_006_G3", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_D12", "KBIMP-_006_G4", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_E12", "KBIMP-_006_G5", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_F12", "KBIMP-_006_G6", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_G12", "KBIMP-_006_G7", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_H10", "KBIMP-_006_G2", SampleID)) 


#### - Filtering seq data for controls and combining different sequence runs ----

#input raw sequence data 

KGLKTK_COI_Domcontigs <- read_tsv(file = "../raw-data2/KGLKTK_2024_OTUDetails.tsv")
CBAYCOIrestplates_Domcontigs <- read_tsv(file = "../raw-data2/CBAY2024_AllPlates_OTUDetails.tsv")
CBAYCOIplate3_Domcontigs <- read_tsv(file = "../raw-data2/Shauna_CBAY2024_Plate3_OTUDetails.tsv")
problemsamples <- read_csv(file = "../raw-data/problemsamples.csv")
View(KGLKTK_COI_Domcontigs)
View(CBAYCOIrestplates_Domcontigs)
View(CBAYCOIplate3_Domcontigs)

length(KGLKTK_COI_Domcontigs$Sample) #168
length(CBAYCOIrestplates_Domcontigs$Sample) #439
length(CBAYCOIplate3_Domcontigs$Sample) #100

#Filtering out/ accoutning for negative controls 
#negative control on plate 3

CBAYCOIplate3_Domcontigs <- CBAYCOIplate3_Domcontigs %>%
  filter(!grepl("10dil", Sample)) %>%
  mutate(ReadCount = ReadCount - ReadCount[Sample == 'NEG_05_H12']) %>%
  filter(!ReadCount <= 0) %>%
  filter(!Sample %in% c("KBIMP_003_G3"))


KGLKTK_COI_allcontigs <- read_tsv(file = "../raw-data/Shauna_COI_KUG2024_TaxonomicAssignments_AllContigs.tsv")

controls_kugluktuk <- KGLKTK_COI_allcontigs %>%
  filter(Sample %in% c("KBIMP-_007_H12", "KBIMP_01_B12")) %>%
  group_by(Plate, Species, Genus) %>%
  summarise(Control_ReadCount = max(ReadCount, na.rm = TRUE), .groups = "drop")

CBAY_COI_allcontigs <- read_tsv(file = "../raw-data/Shauna_COI_CBAY2024_TaxonomicAssignments_AllContigs.tsv")

controls_cbay <- CBAY_COI_allcontigs %>%
  filter(str_detect(Sample, "H12")) %>%
  group_by(Plate, Species, Genus) %>%
  summarise(Control_ReadCount = max(ReadCount, na.rm = TRUE), .groups = "drop")

#negative controls on plate 1,2,4,5,6


CBAYCOIrestplates_Domcontigs <- CBAYCOIrestplates_Domcontigs %>%
  filter(!Sample %in% c("KBIMP_004_H11", "KBIMP_004_D12", "KBIMP_004_E12", "KBIMP-_006_G4",
                        "KBIMP_004_F12", "KBIMP_004_G12", "KBIMP_004_H10", "KBIMP-_006_H1",
                        "KBIMP-_006_H10", "KBIMP-_006_H4", "KBIMP-_006_H5", "KBIMP-_006_H6",
                        "KBIMP-_006_H8", "KBIMP_01_G11", "KBIMP_01_G12", "KBIMP_01_H1", "KBIMP_01_G11", "KBIMP-_006_H11", "KBIMP_002_B11")) %>%
  filter(!Sample %in% problemsamples$`problem samples`) %>%
  left_join(controls_cbay, by = c("Plate", "Species", "Genus")) %>%
  filter(is.na(Control_ReadCount) | ReadCount > Control_ReadCount) %>%
  dplyr::select(-Control_ReadCount) 

#negative controls for Kugluktuk plate 

KGLKTK_COI_Domcontigs <- KGLKTK_COI_Domcontigs %>%
  filter(!Sample %in% c("KBIMP-_008_F4", "KBIMP-_008_F5", "KBIMP_01_B3", "KBIMP_01_F8", "KBIMP_01_B11", "KBIMP_01_D2", "KBIMP_01_E5", "KBIMP_01_F7")) %>%
  left_join(controls_kugluktuk, by = c("Plate", "Species", "Genus")) %>%
  filter(is.na(Control_ReadCount) | ReadCount > Control_ReadCount) %>%
  dplyr::select(-Control_ReadCount) 

length(KGLKTK_COI_Domcontigs$Sample) #156 lost 48 (some of these were due to sampling issues and not controls)
length(CBAYCOIrestplates_Domcontigs$Sample) #359 lost 76 (but this includes the 7 repeats)
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
  dplyr::count(Family, name = "countfamily")

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
  dplyr::count(Genus, name = "countgenera")

ggplot(KBIMP_Domcontigs_bitinggenus, aes(y=countgenera, x=Genus)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 
#2 mosquito genera and 4 black fly genera

rm(KBIMP_Domcontigs_bitinggenus)

#looking at thr black fly species 

KBIMP_Domcontigs_blackflyspecies <- KBIMP_Domcontigs %>%
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

KBIMP_Domcontigs_mosquitoes <- KBIMP_Domcontigs %>%
  filter(Family == "Culicidae") %>%
  select(Species) %>%
  dplyr::count(Species, name = "countspecies")

length(KBIMP_Domcontigs_mosquitoes$Species) #2 species of black fly are seen here but most of them are unknowns 

ggplot(KBIMP_Domcontigs_mosquitoes, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

rm(KBIMP_Domcontigs_mosquitoes, controls_kugluktuk, CBAYCOIplate3_Domcontigs, 
   CBAY_COI_allcontigs, KGLKTK_COI_allcontigs, problemsamples, CBAYCOIrestplates_Domcontigs, controls_cbay)


