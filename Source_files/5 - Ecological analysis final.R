
#To run this data set you need the meatdata from 2024, the fixed site names from 
#2024, the metadata from both places in 2025, the updates species data sets from 
#both years (KBIMP2025_updatedspecies and KBIMP2025_updatedspecies), the data set
#from schafer from 2012 (sr_2012), and the vector dataset from this time 
#this will output multiple plots including the rarefaction curves, vector change plot
#total species richness, and multidimensonal sclaing and well as mulitple statistical tests. 


#loading packages 

library(tidyverse)
library(readr)
library(viridis)
library(parallel) #for running bootstrapping in parallel
library(ggplot2)    #for plotting 
library(betapart)
library(car)
library(DescTools)
library(iNEXT)
library(aRtsy)
library(vegan)
library(ARTool)
library(emmeans)
library(lmerTest)
library(patchwork)
library(brms)
library(ggvenn)

#opening the required data 

CBAY2025_metadata <- read_csv(file = "raw-data2/CBAY2025_metadata.csv")
KGLTK2025_metadata <- read_csv(file = "raw-data2/KGLTK2025_metadata.csv")
condensedsites <- read_csv(file = "raw-data2/condencedsites.csv")
KBIMP2025_updatedspecies <- read_tsv(file = "processed-data/KBIMP2025_updatedspecies.tsv")
KBIMP2024_updatedspecies <- read_tsv(file = "processed-data/KBIMP2024_updatedspecies.tsv")
kbimp2024_sampledata <- read_csv(file = "raw-data2/KBIMP2024_specimendata.csv")
kbimp2024_sitesnamesfixed <- read_csv(file = "raw-data2/KBIMP_meta_sitenamesfixed.csv")
kbimp2024_abundence <- read_csv(file = "raw-data2/KBIMP2024_abundence.csv")
sr_2012 <- read.csv(file = "raw-data2/schafer_2012.csv")
vector_change <- read_csv(file = "raw-data2/vector_change.csv")

#### Invesitgating the number of black flies and mosquitoes from each year ----

#### 2024 ####

#CBAY

sum(kbimp2024_abundence$Sector == "CBAY" ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "CBAY" &
      kbimp2024_abundence$blackfly_abun != 0 & 
      kbimp2024_abundence$mosquito_abun != 0
     ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "CBAY" &
      kbimp2024_abundence$blackfly_abun == 0 & 
      kbimp2024_abundence$mosquito_abun != 0
    ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "CBAY" &
      kbimp2024_abundence$blackfly_abun != 0 & 
      kbimp2024_abundence$mosquito_abun == 0
    ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "CBAY" &
      kbimp2024_abundence$blackfly_abun == 0 & 
      kbimp2024_abundence$mosquito_abun == 0
    ,na.rm = TRUE)

#Kugluktuk

sum(kbimp2024_abundence$Sector == "KGLTK" ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "KGLTK" &
      kbimp2024_abundence$blackfly_abun != 0 & 
      kbimp2024_abundence$mosquito_abun != 0
    ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "KGLTK" &
      kbimp2024_abundence$blackfly_abun == 0 & 
      kbimp2024_abundence$mosquito_abun != 0
    ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "KGLTK" &
      kbimp2024_abundence$blackfly_abun != 0 & 
      kbimp2024_abundence$mosquito_abun == 0
    ,na.rm = TRUE)

sum(kbimp2024_abundence$Sector == "KGLTK" &
      kbimp2024_abundence$blackfly_abun == 0 & 
      kbimp2024_abundence$mosquito_abun == 0
    ,na.rm = TRUE)

##### 2025 #####

CBAY2025_metadata <- CBAY2025_metadata %>%
  filter(!is.na(`Mosquito Head Abundance`))

CBAY2025_metadata$`Blackfly Head Abundance` <- as.numeric(CBAY2025_metadata$`Blackfly Head Abundance`)
CBAY2025_metadata$`Mosquito Head Abundance` <- as.numeric(CBAY2025_metadata$`Mosquito Head Abundance`)
total <- sum(CBAY2025_metadata$`Blackfly Head Abundance`)
total <- sum(CBAY2025_metadata$`Mosquito Head Abundance`)


sum(CBAY2025_metadata$`Blackfly Head Abundance` == 0 &
      CBAY2025_metadata$`Mosquito Head Abundance`!= 0,na.rm = TRUE)

sum(CBAY2025_metadata$`Blackfly Head Abundance` != 0 &
      CBAY2025_metadata$`Mosquito Head Abundance`!= 0,na.rm = TRUE)

sum(CBAY2025_metadata$`Blackfly Head Abundance` == 0 &
      CBAY2025_metadata$`Mosquito Head Abundance`== 0,na.rm = TRUE)

sum(CBAY2025_metadata$`Blackfly Head Abundance` != 0 &
      CBAY2025_metadata$`Mosquito Head Abundance`== 0,na.rm = TRUE)

sum(CBAY2025_metadata$`Blackfly Head Abundance` != 0, na.rm = TRUE)

sum(CBAY2025_metadata$`Mosquito Head Abundance` == 0, na.rm = TRUE)

sum(CBAY2025_metadata$`Mosquito Head Abundance` != 0, na.rm = TRUE)



KGLTK2025_metadata  <- KGLTK2025_metadata %>%
  filter(!is.na(`Mosquito Head Abundance`))

KGLTK2025_metadata $`Blackfly Head Abundance` <- as.numeric(KGLTK2025_metadata $`Blackfly Head Abundance`)
KGLTK2025_metadata $`Mosquito Head Abundance` <- as.numeric(KGLTK2025_metadata $`Mosquito Head Abundance`)
total <- sum(KGLTK2025_metadata$`Blackfly Head Abundance`)
total <- sum(KGLTK2025_metadata$`Mosquito Head Abundance`)


sum(KGLTK2025_metadata$`Mosquito Head Abundance` != 0 &
      KGLTK2025_metadata$`Blackfly Head Abundance`!= 0 , na.rm = TRUE)

sum(KGLTK2025_metadata$`Mosquito Head Abundance` != 0 &
      KGLTK2025_metadata$`Blackfly Head Abundance`== 0 , na.rm = TRUE)

sum(KGLTK2025_metadata$`Mosquito Head Abundance` == 0 &
      KGLTK2025_metadata$`Blackfly Head Abundance`!= 0 , na.rm = TRUE)

sum(KGLTK2025_metadata$`Mosquito Head Abundance` == 0 &
      KGLTK2025_metadata$`Blackfly Head Abundance`== 0 , na.rm = TRUE)

sum(KGLTK2025_metadata$`Blackfly Head Abundance` == 0, na.rm = TRUE)
sum(KGLTK2025_metadata$`Blackfly Head Abundance` != 0, na.rm = TRUE)
sum(KGLTK2025_metadata$`Mosquito Head Abundance` == 0, na.rm = TRUE)
sum(KGLTK2025_metadata$`Mosquito Head Abundance` != 0, na.rm = TRUE)


#### - combining 2024 and 2025 data into one data set ----

#fixing sample names and selecting for required columns 

KBIMP2025_updatedclean <- KBIMP2025_updatedspecies %>%
  filter(Sample != "Outgroup", Family %in% c("Culicidae", "Simuliidae")) %>%
  mutate(Sample = str_extract(Sample, "^[A-Za-z]+_?\\d+"),
         Sample = str_replace(Sample, "_", "")) %>%
  select(Sample, Species, Family) 

KBIMP2024_updatedclean <- KBIMP2024_updatedspecies %>%
  full_join(kbimp2024_sampledata, join_by(Sample == SampleID)) %>% 
  filter(Sample != "Outgroup", Family %in% c("Culicidae", "Simuliidae")) %>%
  select(FieldID, Species, Family) %>%
  dplyr::rename(Sample = FieldID) 

#combining by stacking rows 

KBIMP_combined <- KBIMP2025_updatedclean %>%
  bind_rows(KBIMP2024_updatedclean) %>%
  distinct()

#### - Preparing the site metadata for analysis ----

##### 2024 sample data #####

#looking at the site names

KBIMP_2024_siteslook <- kbimp2024_sampledata %>%
  group_by(Site) %>%
  summarise(countsites = n_distinct(FieldID)) 

ggplot(KBIMP_2024_siteslook, aes(y= countsites, x= Site)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(KBIMP_2024_siteslook)

#combine good sites sheet with other sampledata sheet and adding in the sample month

kbimp2024_sampledata <- kbimp2024_sampledata %>%
  inner_join(kbimp2024_sitesnamesfixed, join_by(FieldID == FieldID)) %>%
  mutate(
    date_parts = str_split(`Collection Date`, "-"),
    end_date_raw = sapply(date_parts, `[`, 2),
    end_date_raw = ifelse(is.na(end_date_raw), `Collection Date`, end_date_raw),
    Date_parsed = parse_date_time(end_date_raw,
                                  orders = c("ymd", "mdy", "dmy", "dm", "md", "m")),
    Month = month(Date_parsed, label = TRUE))

#looking at the new site names 

KBIMP_2024_siteslook <- kbimp2024_sampledata %>%
  group_by(ExactSite) %>%
  summarise(countsites = n_distinct(FieldID)) 


ggplot(KBIMP_2024_siteslook, aes(y= countsites, x= ExactSite)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(KBIMP_2024_siteslook)

#fixing the samples which were switched during extractions (based on lab notes)

kbimp2024_sampledata <- kbimp2024_sampledata %>%
  mutate(SampleID= gsub("KBIMP_004_H11", "KBIMP-_006_G3", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_D12", "KBIMP-_006_G4", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_E12", "KBIMP-_006_G5", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_F12", "KBIMP-_006_G6", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_G12", "KBIMP-_006_G7", SampleID)) %>%
  mutate(SampleID= gsub("KBIMP_004_H10", "KBIMP-_006_G2", SampleID)) 

##### prearing the 2025 data #####

CBAY2025_metadata <- CBAY2025_metadata %>%
  mutate(Date_parsed = parse_date_time(`Date collected`,
                                       orders = c("ymd", "mdy", "dmy")),
         Month = month(Date_parsed, label = TRUE))%>%
  select(Month, Sample, Site, `Sample Replicate`, `Sample type collection method`,
         Lat, Lon, `Habitat type`, `Distance from water (m)`, `Mosquito Head Abundance`, 
         `Blackfly Head Abundance`, `Temp C (Sweep)`, `Relative Humidity (Sweep)`) %>%
  right_join(condensedsites)

KGLTK2025_metadata <- KGLTK2025_metadata %>%
  mutate(Date_parsed = parse_date_time(`Date set`,
                                       orders = c("ymd", "mdy", "dmy", "dm")),
         Month = month(Date_parsed, label = TRUE))

#setting up metadata for analysis 

KGLTK2025_metadata <- KGLTK2025_metadata %>%
  select(Sample, `Sample type collection method`, Month) 

KBIMP2025_metadata <- CBAY2025_metadata %>% 
  select(Sample, `Sample type collection method`, Month) %>%
  rbind(KGLTK2025_metadata) %>%
  filter(!is.na(Sample)) %>%
  mutate(Sector = str_extract(Sample, "^[A-Za-z]+")) %>%
  mutate(Year = 2025) %>%
  dplyr::rename(SamplingMethod = `Sample type collection method`) 

##### combining 2024 and 2025 data into one metadata file ----- 

meta_data <- kbimp2024_sampledata %>%
  dplyr::rename(Sample = FieldID) %>%
  dplyr::rename(Lat = Lat.x) %>%
  select(Sample, SamplingMethod, Sector, Month) %>%
  mutate(Year = 2024) %>%
  bind_rows(KBIMP2025_metadata) %>%
  distinct() %>%
  mutate(SamplingMethod = case_when(
    grepl("malaise", SamplingMethod, ignore.case = TRUE) ~ "Malaise Trap",
    grepl("sweep", SamplingMethod, ignore.case = TRUE) ~ "Sweep Net",
    grepl("aspirator|people", SamplingMethod, ignore.case = TRUE) ~ "Aspirator",
    TRUE ~ SamplingMethod)) 

sum(meta_data$Month == "Jul" &
      meta_data$Sector == "CBAY",
    na.rm = TRUE)

sum(meta_data$Month == "Jul" &
      meta_data$Sector == "KGLTK",
    na.rm = TRUE)

sum(meta_data$Month == "Aug" &
      meta_data$Sector == "CBAY",
    na.rm = TRUE)

sum(meta_data$Month == "Aug" &
      meta_data$Sector == "KGLTK",
    na.rm = TRUE)

####  Making iNEXT graph ----

##### iNEXT both places, both families #####  

iNEXT <- KBIMP_combined %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 233,
          KGLTK = 122, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext <- iNEXT(iNEXT, q=0, datatype="incidence_freq")

em.inext$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_total <- as.data.frame(em.inext$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) 

em.inext$DataInfo

##### iNEXT both places, just black flies both months #####  

iNEXT_sim <- KBIMP_combined %>%
  filter(Family == "Simuliidae" )%>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 233,
          KGLTK = 122, .before = 1) %>%
  column_to_rownames(var = "Species")


em.inext.sim <- iNEXT(iNEXT_sim, q=0, datatype="incidence_freq")

em.inext.sim$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_bf <- as.data.frame(em.inext.sim$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(Family = "Simuliidae") %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) 


em.inext.sim$DataInfo

iNext_bf <- ggiNEXT(em.inext.sim)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18))

ggsave("plots/inext2025bf.png", iNext_bf , width = 6, height = 5, dpi = 300, bg = "transparent")

##### Both places, mosquitoes, both months #####

iNEXT_cul <- KBIMP_combined %>%
  filter(Family == "Culicidae")%>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 233,
          KGLTK = 122, .before = 1) %>%
  column_to_rownames(var = "Species")


em.inext.cul <- iNEXT(iNEXT_cul, q=0, datatype="incidence_freq")

em.inext.cul$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_mos <- as.data.frame(em.inext.cul$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(Family = "Culicidae") %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) 

em.inext.cul$DataInfo

iNext_mos <- ggiNEXT(em.inext.cul)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("#000099", "#FFC000")) +
  scale_fill_manual( values = c("#000099", "#FFC000")) +
  scale_shape_manual( values = c(15, 16, 17, 18))

ggsave("plots/inext2025mos.png", iNext_mos , width = 6, height = 4, dpi = 300, bg = "transparent")

##### Both places, just mosquitoes, just july #####

Julsamples <- meta_data %>%
  filter(Month == "Jul")

iNEXT_cul_jul <- KBIMP_combined %>%
  filter(Family == "Culicidae")%>%
  filter(Sample %in% Julsamples$Sample) %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 127,
          KGLTK = 32, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext.cul.jul <- iNEXT(iNEXT_cul_jul, q=0, datatype="incidence_freq")

em.inext.cul.jul$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_mos_jul <- as.data.frame(em.inext.cul.jul$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(Family = "Culicidae") %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(Month = "July")

##### Both places, just blackflies, just july #####

iNEXT_bf_jul <- KBIMP_combined %>%
  filter(Family == "Simuliidae")%>%
  filter(Sample %in% Julsamples$Sample) %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 127,
          KGLTK = 32, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext.bf.jul <- iNEXT(iNEXT_bf_jul, q=0, datatype="incidence_freq")

em.inext.bf.jul$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_bf_jul <- as.data.frame(em.inext.bf.jul$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(Family = "Simuliidae") %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(Month = "July")

em.inext.bf.jul$DataInfo

##### Both places, just mosquitoes, just august #####

Augsamples <- meta_data %>%
  filter(Month == "Aug")

iNEXT_cul_aug <- KBIMP_combined %>%
  filter(Family == "Culicidae")%>%
  filter(Sample %in% Augsamples$Sample) %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 84,
          KGLTK = 54, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext.cul.aug <- iNEXT(iNEXT_cul_aug, q=0, datatype="incidence_freq")

em.inext.cul$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_mos_aug <- as.data.frame(em.inext.cul.aug$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(Family = "Culicidae") %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(Month = "August")

em.inext.cul$DataInfo


##### Both places, just blackflies, just august #####

iNEXT_bf_aug <- KBIMP_combined %>%
  filter(Family == "Simuliidae")%>%
  filter(Sample %in% Augsamples$Sample) %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 84,
          KGLTK = 54, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext.bf.aug <- iNEXT(iNEXT_bf_aug, q=0, datatype="incidence_freq")

em.inext.bf.aug$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_bf_aug <- as.data.frame(em.inext.bf.aug$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(Family = "Simuliidae") %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(Month = "August")

em.inext.cul$DataInfo

#### - Investigating differences since 2012 ----

##### venn diagram/ determining which species are different #####

sr_2012 <- sr_2012 %>%
  mutate(taxon = if_else(taxon == "Aedes nigripes", "Aedes nigripes/Aedes impiger", taxon)) %>%
  filter(!taxon == "Aedes impiger")

venndiagramspecieslist <- KBIMP_combined  %>%
  mutate(Species = str_replace(Species, "Aedes punctor/Aedes hexodontus", "Aedes hexodontus")) %>%
  mutate(Species = if_else(Species == "Aedes nigripes/impiger", "Aedes nigripes/Aedes impiger", Species)) %>%
  separate_rows(Species, sep = "/") %>%
  mutate(Species = str_replace(Species, "Simulium arcticum complex sp", "Simulium arcticum complex")) %>%
  mutate(Species = str_replace(Species, "Simulium verecundum complex sp", "Simulium verecundum complex")) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  mutate(presence = 1) %>%
  select(-Sample, -Family) %>%
  distinct() %>%
  group_by(Species, region) %>% 
  summarise(across(everything(), ~ as.integer(any(.x == 1))), .groups = "drop") %>%
  ungroup() %>%
  pivot_wider(names_from = "region", values_from = "presence") %>%
  full_join(sr_2012, join_by(Species == taxon)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) 

venncbay <- venndiagramspecieslist %>%
  select(CBAY, CbB, Species) %>%
  column_to_rownames(var = "Species") %>%
  filter(rowSums(across(everything())) > 0) %>%
  rownames_to_column(var = "Species") %>%
  pivot_longer(-Species, names_to = "sector", values_to = "presence") %>%
  filter(presence == 1) %>%
  select(sector, Species) %>%
  group_by(sector) %>%
  summarise(species_list = list(unique(Species))) %>%
  deframe()

setdiff(venncbay$CBAY, venncbay$CbB)
setdiff(venncbay$CbB, venncbay$CBAY)

names(venncbay) <- c("CBAY" = "Cambridge Bay 2024and 2025", "CbB" = "Cambridge Bay 2012")

vennkug <- venndiagramspecieslist %>%
  select(KGLTK,Kug, Species) %>%
  column_to_rownames(var = "Species") %>%
  filter(rowSums(across(everything())) > 0) %>%
  rownames_to_column(var = "Species") %>%
  pivot_longer(-Species, names_to = "sector", values_to = "presence") %>%
  filter(presence == 1) %>%
  select(sector, Species) %>%
  group_by(sector) %>%
  summarise(species_list = list(unique(Species))) %>%
  deframe()

setdiff(vennkug$KGLTK, vennkug$Kug)
setdiff(vennkug$Kug, vennkug$KGLTK)

ggsave("plots/vennkug.png", vennkug, width = 4, height = 2, dpi = 300)

##### making the figure for the change in the number of vectors #####

#getting the confidence intervals for vectors 

iNEXT_vectors <- KBIMP_combined %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  filter(Species %in% vector_change$Species) %>%
  add_row(Species = "sampling_extent",
          CBAY = 233,
          KGLTK = 122, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext_vectors <- iNEXT(iNEXT_vectors, q=0, datatype="incidence_freq")

em.inext_vectors$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_vectors <- as.data.frame(em.inext_vectors$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) 

em.inext_vectors$DataInfo

#getting the confidence intervals for non-vectors 

iNEXT_nonvectors <- KBIMP_combined %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  filter(!Species %in% vector_change$Species) %>%
  add_row(Species = "sampling_extent",
          CBAY = 233,
          KGLTK = 122, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext_nonvectors <- iNEXT(iNEXT_nonvectors, q=0, datatype="incidence_freq")

em.inext_nonvectors$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_nonvectors <- as.data.frame(em.inext_nonvectors$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) 

em.inext_nonvectors$DataInfo

#getting the vector, non-vector, total data from the 2012 dataset

  
vector2012 <- sr_2012 %>%
  
  pivot_longer(
    cols = -taxon,
    names_to = "Site",
    values_to = "Presence") %>%
  
  filter(Presence == 1) %>%

  mutate(InList = taxon %in% vector_change$Species) %>%
  
  group_by(Site) %>%
  summarise(
    Total_species = n(),
    vector_species = sum(InList),
    nonvector_species = sum(!InList))


#figure with both years for cbay

vectorchange <- ggplot(vector_change, aes(x = Year, y =  Count, group = Type, colour =  Type)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(
    aesthetics = c("colour", "fill"),
    values = c("Vector" = "#000099", "Non-Vector" =  "#FFC000", "Total" = "Black"),
    name = "Type") +
  theme_bw(base_size = 15) +
  xlab("Year") +
  ylab("Total Species") +
  theme(
    axis.text.x = element_text(angle = 0, size = 10, color = "black"),
    axis.text.y = element_text(angle = 0, size = 10, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.background = element_rect(fill = NA, color = NA), legend.background = element_rect(fill = NA, color = NA), 
    strip.background = element_rect(fill = NA, color = NA),
    strip.text  = element_text(face = "bold", size = 14)) 


ggsave("plots/vectorchangecbay.png", vectorchange , width = 8, height = 4, dpi = 300, bg = "transparent")

#figure with both years for kugluktuk 

vectorchangekug <- ggplot(vector_change_kug, aes(x = Year, y =  Count, group = Type, colour =  Type)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(
    aesthetics = c("colour", "fill"),
    values = c("Vector" = "#000099", "Non-Vector" =  "#FFC000", "Total" = "Black"),
    name = "Type") +
  theme_bw(base_size = 15) +
  xlab("Year") +
  ylab("Total Species") +
  theme(
    axis.text.x = element_text(angle = 0, size = 10, color = "black"),
    axis.text.y = element_text(angle = 0, size = 10, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.background = element_rect(fill = NA, color = NA), legend.background = element_rect(fill = NA, color = NA), 
    strip.background = element_rect(fill = NA, color = NA),
    strip.text  = element_text(face = "bold", size = 14)) 


ggsave("plots/vectorchangekug.png", vectorchangekug , width = 8, height = 4, dpi = 300, bg = "transparent")


#### - Alpha div analysis ----


#create data set by counting the number of each species in each sample

alpha_en_data <- meta_data %>%
  filter(SamplingMethod %in% c("Malaise Trap", "Sweep Net")) 

families <- c("Culicidae", "Simuliidae") 

speciesrich_CBAYvs_KGLTK <- KBIMP_combined %>%
  filter(Sample %in% alpha_en_data$Sample) %>%
  inner_join(alpha_en_data) %>%
  filter(!Month %in% c("Jun", "Sep")) %>%
  group_by(Sector, Sample, Month, Year, Family) %>%
  summarise(Speciessum = n_distinct(Species, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>%
  group_by(Sector, Sample, Year, Month) %>%  
  complete(
    Family = families,
    fill = list(Speciessum = 0)
  ) %>%
  ungroup() %>%
  filter(!Sample == "KGLTK") %>%
  filter(!is.na(Family))

#exploring data for normality and equality of variance

leveneTest(Speciessum ~ Sector, data = speciesrich_CBAYvs_KGLTK) #p-value greater than 0.05 - vairances equal 
shapiro.test(speciesrich_CBAYvs_KGLTK$Speciessum) #p-value less than 0.05 - data not normal

#running a non-parmetric test because the data did not meet the above requirements

m <- art(Speciessum ~ factor(Sector) * factor(Month) * factor(Family) + Error(factor(Year)),
         data = speciesrich_CBAYvs_KGLTK)
anova(m)

modeluniqspsample <- kruskal.test(Speciessum ~ Sector, data = speciesrich_CBAYvs_KGLTK)

modeluniqspsample


#creating data set of total per month for plot

confidence_month <- bind_rows(confidenceinterval_mos_jul, 
                              confidenceinterval_bf_jul, 
                              confidenceinterval_mos_aug,
                              confidenceinterval_bf_aug)

confidence_month$Month <- factor(
  confidence_month$Month,
  levels = c("July", "August")) 

#creating data set with total over all and confidence intverals for plot

confidence <- bind_rows(confidenceinterval_bf, confidenceinterval_mos)

#creating the plot

speciesrichplot <- ggplot() +
  
  geom_line(data = confidence_month,
    aes(x = Month, y = qD, colour = Sector, group = Sector)) +
  
  geom_point(data = confidence_month, aes(x = Month, y = qD, colour = Sector),
    size = 3) +
  
  geom_errorbar(data = confidence_month, aes(x = Month, ymin = qD.LCL, ymax = qD.UCL, colour = Sector),
                width = 0.15) +
  
  geom_point(data = confidence, aes(x = x, y = qD, colour = Sector),
    shape = 17, size = 3) +
  
  geom_errorbar(data = confidence, aes(x = x, ymin = qD.LCL, ymax = qD.UCL, colour = Sector),
    width = 0.15) +
  
  scale_colour_manual(
    values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"),
    breaks = c("CBAY", "KGLTK"),
    labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n",
                "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"), name = "Region") +
  
  theme_bw(base_size = 15) +
  xlab("Sampling Month") +
  ylab("Species Richness") +
  
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black"),
    axis.text.y = element_text(angle = 0, size = 10, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.background = element_rect(fill = NA, color = NA), legend.background = element_rect(fill = NA, color = NA), 
    strip.background = element_rect(fill = NA, color = NA),
    strip.text  = element_text(face = "bold", size = 14)) +
  
  facet_wrap(~Family) 

ggsave("plots/lineSRtotalSRcbayvskug.png", speciesrichplot , width = 10, height = 4, dpi = 300, bg = "transparent")

#### - Multidimentional 2024 and 2025 ----

##### Just mosquitoes -----

KBIMP_speciesmatrix_mosquitoes <- KBIMP_combined %>%
  filter(Family == "Culicidae") %>%
  select(Species, Sample) %>%
  dplyr::count(Sample, Species, name = "Abundance") %>% 
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

multi_en_data_mos <- meta_data %>%
  filter(Sample %in% rownames(KBIMP_speciesmatrix_mosquitoes)) %>%
  distinct() %>%
  filter(SamplingMethod %in% c("Malaise Trap","Sweep Net")) 

#making sure only the samples we filtered above are in the species matrix 

KBIMP_speciesmatrix_mosquitoes <- KBIMP_speciesmatrix_mosquitoes %>%
  filter(rownames(KBIMP_speciesmatrix_mosquitoes) %in% multi_en_data_mos$Sample) 

# Order species matrix

KBIMP_speciesmatrix_mosquitoes <- KBIMP_speciesmatrix_mosquitoes[order(rownames(KBIMP_speciesmatrix_mosquitoes)), ]

# Order environmental data

multi_en_data_mos <- multi_en_data_mos %>% arrange(Sample) 

# Bray-Curtis dissimilarity

bray_dist <- vegdist(KBIMP_speciesmatrix_mosquitoes, method = "bray")

#permanova 

adonis2(bray_dist ~ Sector + Year + SamplingMethod + Month,
        data = multi_en_data_mos, permutations = 999, by="margin")

#multidimensonal scaling

NMDS <- metaMDS(KBIMP_speciesmatrix_mosquitoes, k=2)

NMDS$stress

#anosim for Sector

anosim(KBIMP_speciesmatrix_mosquitoes, multi_en_data_mos$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Year

anosim(KBIMP_speciesmatrix_mosquitoes, multi_en_data_mos$Year ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Sampling method 

anosim(KBIMP_speciesmatrix_mosquitoes, multi_en_data_mos$SamplingMethod ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Month

anosim(KBIMP_speciesmatrix_mosquitoes, multi_en_data_mos$Month ,
       permutations = 999,distance = "bray", strata = NULL)

# Fit environmental vector for CollectionWeek
env_fit <- envfit(NMDS, multi_en_data_mos, permutations = 999)


#extracting the nmds data into a dataframe 

nmds_scores <- as.data.frame(vegan::scores(NMDS, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  full_join(multi_en_data_mos)
en_coord_cat <- as.data.frame(vegan::scores(env_fit, "factors")) * ordiArrowMul(env_fit)

# Apply the function to each group to create the polygons on the ggplot nmds 
hulls <- nmds_scores %>%
  group_by(Sector) %>%
  group_split() %>%
  lapply(calculate_hull) %>%
  bind_rows()

#making the ggplot nmds plot

nmdsplotmos <- ggplot(data = nmds_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_polygon(data = hulls,aes(x = NMDS1, y = NMDS2, group = Sector, fill = Sector),
               alpha = 0.3, colour = "black") +
  
  geom_point(data = nmds_scores, aes(colour = Sector), size = 2, alpha = 1) +
  
  scale_colour_manual(aesthetics = c("colour", "fill"), values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"),
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 9, colour = "black"), legend.position = "top")


nmdsplotmos 


ggsave("plots/nmdsplotmos2024and2025.png", nmdsplotmos, width = 3, height = 3, dpi = 300)


##### just black flies -----

KBIMP_speciesmatrix_blackflies <- KBIMP_combined %>%
  filter(!Family == "Culicidae") %>%
  select(Species, Sample) %>%
  dplyr::count(Sample, Species, name = "Abundance") %>% 
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

multi_en_data_bf <- meta_data %>%
  filter(Sample %in% rownames(KBIMP_speciesmatrix_blackflies)) %>%
  distinct() %>%
  filter(SamplingMethod %in% c("Malaise Trap","Sweep Net")) 

KBIMP_speciesmatrix_blackflies <- KBIMP_speciesmatrix_blackflies %>%
  filter(rownames(KBIMP_speciesmatrix_blackflies) %in% multi_en_data_bf$Sample) 

# Order species matrix

KBIMP_speciesmatrix_blackflies <- KBIMP_speciesmatrix_blackflies[order(rownames(KBIMP_speciesmatrix_blackflies)), ]

# Order environmental data

multi_en_data_bf <- multi_en_data_bf %>% arrange(Sample) 

# Bray-Curtis dissimilarity

bray_dist <- vegdist(KBIMP_speciesmatrix_blackflies, method = "bray")

adonis2(bray_dist ~ Sector + Year + SamplingMethod + Month,
        data = multi_en_data_bf, permutations = 999, by="margin")

NMDS <- metaMDS(KBIMP_speciesmatrix_blackflies, k=3)

NMDS$stress

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$Year ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$SamplingMethod ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$Month ,
       permutations = 999,distance = "bray", strata = NULL)

# Fit environmental vector for CollectionWeek
env_fit <- envfit(NMDS, multi_en_data_bf, permutations = 999)

#convert the data from nmds into a dataframe for plotting 

nmds_scores <- as.data.frame(vegan::scores(NMDS, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  full_join(multi_en_data_bf)
en_coord_cat <- as.data.frame(vegan::scores(env_fit, "factors")) * ordiArrowMul(env_fit)


# Apply the function to each group to create the polygons on the ggplot nmds 
hulls <- nmds_scores %>%
  group_by(Sector) %>%
  group_split() %>%
  lapply(calculate_hull) %>%
  bind_rows()

#making the ggplot nmds plot

nmdsplotbf <- ggplot(data = nmds_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, group = Sector, fill = Sector),
               alpha = 0.3, colour = "black") +
  
  geom_point(data = nmds_scores, aes(colour = Sector), size = 2, alpha = 1) +
  
  geom_text(aes(label = Sample), vjust = -0.5) + 
  
  scale_colour_manual(aesthetics = c("colour", "fill"), values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 9, colour = "black"), legend.position = "top")


nmdsplotbf


ggsave("plots/nmdsplotbf2024and2025.png", nmdsplotbf, width = 3, height = 3, dpi = 300)

#two of the points are making it really hard to see the data so I am going to remove these points for plotting 
##### - just black flies with those two kugluktuk samples removed -----

KBIMP_speciesmatrix_blackflies2 <- KBIMP_combined %>%
  filter(!Family == "Culicidae") %>%
  select(Species, Sample) %>%
  dplyr::count(Sample, Species, name = "Abundance") %>% 
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

multi_en_data_bf <- meta_data %>%
  filter(Sample %in% rownames(KBIMP_speciesmatrix_blackflies)) %>%
  distinct() %>%
  filter(SamplingMethod %in% c("Malaise Trap","Sweep Net")) %>%
  filter(!Sample %in% c("KGLTK0137", "KGLTK0103"))

KBIMP_speciesmatrix_blackflies2 <- KBIMP_speciesmatrix_blackflies2 %>%
  filter(rownames(KBIMP_speciesmatrix_blackflies2) %in% multi_en_data_bf$Sample) 

# Order species matrix
KBIMP_speciesmatrix_blackflies2 <- KBIMP_speciesmatrix_blackflies2[order(rownames(KBIMP_speciesmatrix_blackflies2)), ]

# Order environmental data
multi_en_data_bf <- multi_en_data_bf %>% arrange(Sample) 

# Bray-Curtis dissimilarity
bray_dist <- vegdist(KBIMP_speciesmatrix_blackflies2, method = "bray")

adonis2(bray_dist ~ Sector + Year + SamplingMethod + Month,
        data = multi_en_data_bf, permutations = 999, by="margin")

NMDS <- metaMDS(KBIMP_speciesmatrix_blackflies2, k=3)

NMDS$stress

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$Year ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$SamplingMethod ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$Month ,
       permutations = 999,distance = "bray", strata = NULL)

# Fit environmental vector for CollectionWeek
env_fit <- envfit(NMDS, multi_en_data_bf, permutations = 999)

#converting the nmds data to a dataframe for plotting
nmds_scores <- as.data.frame(vegan::scores(NMDS, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  full_join(multi_en_data_bf)
en_coord_cat <- as.data.frame(vegan::scores(env_fit, "factors")) * ordiArrowMul(env_fit)

# Apply the function to each group to create the polygons on the ggplot nmds 
hulls <- nmds_scores %>%
  group_by(Sector) %>%
  group_split() %>%
  lapply(calculate_hull) %>%
  bind_rows()

#making the ggplot nmds plot

nmdsplotbf2 <- ggplot(data = nmds_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, group = Sector, fill = Sector),
               alpha = 0.3, colour = "black") +
  
  geom_point(data = nmds_scores, aes(colour = Sector), size = 2, alpha = 1) +
  
  scale_colour_manual(aesthetics = c("colour", "fill"), values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"),
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n",
                                                               "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"), name = "Region") +
  
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 9, colour = "black"), legend.position = "top")


nmdsplotbf2 

ggsave("plots/nmdsplotbf22024and2025.png", nmdsplotbf2, width = 3, height = 3, dpi = 300)

##### mosquitoes and black flies just sweep nets -----

KBIMP_speciesmatrix_everything <- KBIMP_combined %>%
  select(Species, Sample) %>%
  dplyr::count(Sample, Species, name = "Abundance") %>% 
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

multi_en_data_swn <- meta_data %>%
  filter(Sample %in% rownames(KBIMP_speciesmatrix_everything)) %>%
  distinct() %>%
  filter(SamplingMethod %in% c("Sweep Net")) 

KBIMP_speciesmatrix_everything <- KBIMP_speciesmatrix_everything %>%
  filter(rownames(KBIMP_speciesmatrix_everything) %in% multi_en_data_swn$Sample) 

# Order species matrix
KBIMP_speciesmatrix_everything <- KBIMP_speciesmatrix_everything[order(rownames(KBIMP_speciesmatrix_everything)), ]

# Order environmental data
multi_en_data_swn <- multi_en_data_swn %>% arrange(Sample) 

# Bray-Curtis dissimilarity
bray_dist <- vegdist(KBIMP_speciesmatrix_everything, method = "bray")

adonis2(bray_dist ~ Sector + Year + Month,
        data = multi_en_data_swn, permutations = 999, by="margin")

NMDS <- metaMDS(KBIMP_speciesmatrix_everything, k=4)

NMDS$stress

anosim(KBIMP_speciesmatrix_everything, multi_en_data_swn$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_everything, multi_en_data_swn$Year ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_everything, multi_en_data_swn$Month ,
       permutations = 999,distance = "bray", strata = NULL)

##### mosquitoes and black flies just malaise traps  -----

KBIMP_speciesmatrix_everything <- KBIMP_combined %>%
  select(Species, Sample) %>%
  dplyr::count(Sample, Species, name = "Abundance") %>% 
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

multi_en_data_mat <- meta_data %>%
  filter(Sample %in% rownames(KBIMP_speciesmatrix_everything)) %>%
  distinct() %>%
  filter(SamplingMethod %in% c("Malaise Trap")) 

KBIMP_speciesmatrix_everything <- KBIMP_speciesmatrix_everything %>%
  filter(rownames(KBIMP_speciesmatrix_everything) %in% multi_en_data_mat$Sample) 

# Order species matrix
KBIMP_speciesmatrix_everything <- KBIMP_speciesmatrix_everything[order(rownames(KBIMP_speciesmatrix_everything)), ]

# Order environmental data
multi_en_data_mat <- multi_en_data_mat %>% arrange(Sample) 

# Bray-Curtis dissimilarity
bray_dist <- vegdist(KBIMP_speciesmatrix_everything, method = "bray")

adonis2(bray_dist ~ Sector + Year + Month,
        data = multi_en_data_mat, permutations = 999, by="margin")

NMDS <- metaMDS(KBIMP_speciesmatrix_everything, k=5)

NMDS$stress

anosim(KBIMP_speciesmatrix_everything, multi_en_data_mat$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_everything, multi_en_data_mat$Year ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_everything, multi_en_data_mat$Month ,
       permutations = 999,distance = "bray", strata = NULL)

##### both sampling methods mosqutioes and black flies combined  -----

KBIMP_speciesmatrix_everything <- KBIMP_combined %>%
  select(Species, Sample) %>%
  dplyr::count(Sample, Species, name = "Abundance") %>% 
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

multi_en_data <- meta_data %>%
  filter(Sample %in% rownames(KBIMP_speciesmatrix_everything)) %>%
  distinct() %>%
  filter(SamplingMethod %in% c("Malaise Trap", "Sweep Net")) 

KBIMP_speciesmatrix_everything <- KBIMP_speciesmatrix_everything %>%
  filter(rownames(KBIMP_speciesmatrix_everything) %in% multi_en_data$Sample) 

# Order species matrix
KBIMP_speciesmatrix_everything <- KBIMP_speciesmatrix_everything[order(rownames(KBIMP_speciesmatrix_everything)), ]

# Order environmental data
multi_en_data <- multi_en_data %>% arrange(Sample) 

# Bray-Curtis dissimilarity
bray_dist <- vegdist(KBIMP_speciesmatrix_everything, method = "bray")

adonis2(bray_dist ~ Sector + Year + Month + SamplingMethod,
        data = multi_en_data, permutations = 999, by="margin")

NMDS <- metaMDS(KBIMP_speciesmatrix_everything, k=4)

NMDS$stress

anosim(KBIMP_speciesmatrix_everything, multi_en_data$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_everything, multi_en_data$Year ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_everything, multi_en_data$SamplingMethod ,
       permutations = 999,distance = "bray", strata = NULL)

anosim(KBIMP_speciesmatrix_everything, multi_en_data$Month ,
       permutations = 999,distance = "bray", strata = NULL)

#### - PART Beta diveristy comparing the two communties ----

species_matrix_communties <- KBIMP_speciesmatrix_everything %>%
  rownames_to_column(var = "Sample") %>%
  full_join(multi_en_data, join_by(Sample == Sample)) %>%
  select(-Year, -SamplingMethod, -Month) %>%
  group_by(Sector) %>%
  summarise(across(-Sample, max), .groups = "drop") %>%
  column_to_rownames(var = "Sector")

#run beatdiverity analysis 

betacommunties <- beta.pair(species_matrix_communties, index.family = "sorensen")  

print(betacommunties$beta.sim)
print(betacommunties$beta.sne)
print(betacommunties$beta.sor)

#### - Number of insects etc ---- 




