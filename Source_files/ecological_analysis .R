#### - Cleaning the data ----


#load required packages

library(readr)
library(tidyverse)
library(ggplot2)
library(car)
library(vegan)
library(betapart)
library(reshape2)
library(DescTools)

#### PART 1 - Cleaning the COI  data ----

#open data of the dominant contigs into R

CBAYCOIrestplates_Domcontigs <- read_tsv(file = "../raw-data/Shauna_COI_CBAY2024_TaxonomicAssignments_DominantContigs.tsv")
CBAYCOIplate3_Domcontigs <- read_tsv(file = "../raw-data/Shauna_100dilution_TaxonomicAssignments_DominantContigs.tsv")

View(CBAYCOIrestplates_Domcontigs)
View(CBAYCOIplate3_Domcontigs)

#combining plate 3 with the rest of the plates

CBAYCOI_Domcontigs <- rbind(CBAYCOIrestplates_Domcontigs, CBAYCOIplate3_Domcontigs)

View(CBAYCOI_Domcontigs)

#lets look at the amount of wrong IDs (or IDs that aren't black flies/mosquitoes)

#looking at it by Family
CBAYCOI_Domcontigs_family <- CBAYCOI_Domcontigs %>%
  select(Family) %>%
  count(Family, name = "countfamily")
view(CBAYCOI_Domcontigs_family)

ggplot(CBAYCOI_Domcontigs_family, aes(y=countfamily, x=Family)) +
  geom_col(fill = "skyblue3") +
  theme_bw()  #7 families which weren't mosquitoes or black flies 

rm(CBAYCOI_Domcontigs_family)

#now lets look at the genus diversity of only black flies and mosquitoes

CBAYCOI_Domcontigs_bitinggenus <- CBAYCOI_Domcontigs %>%
  filter(Family != "Empididae") %>%
  filter(Family != "Muscidae") %>%
  filter(Family != "Chironomidae") %>%
  filter(Family != "Dolichopodidae") %>%
  filter(Family != "Scathophagidae") %>%
  filter(Family != "unknown") %>%
  select(Genus) %>%
  count(Genus, name = "countgenera")

ggplot(CBAYCOI_Domcontigs_bitinggenus, aes(y=countgenera, x=Genus)) +
  geom_col(fill = "skyblue3") +
  theme_bw() 
#all mosqutioes were aedes - but we saw a split in black flies between 4 different genera

rm(CBAYCOI_Domcontigs_bitinggenus)

#it seems like all the mosquitoes were only identified to the genus level so lets look at the species of black flies

CBAYCOI_Domcontigs_blackflyspecies <- CBAYCOI_Domcontigs %>%
  filter(Family == "Simuliidae") %>%
  select(Species) %>%
  count(Species, name = "countspecies")

ggplot(CBAYCOI_Domcontigs_blackflyspecies, aes(y=countspecies, x=Species)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we got 7 different species plus unknown species

rm(CBAYCOI_Domcontigs_blackflyspecies)

#### PART 2 - Combing with sample metadata for analysis of sampling sites ----

CBAY2024_sampledata <- read_csv(file = "../raw-data/CBAY2024_sampledata.csv")

view(CBAY2024_sampledata)

#Combing different sites based on variation in sample names and coordinate location

CBAY2024_sites <- CBAY2024_sampledata %>%
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
  mutate(ExactSite= gsub("creek by long pt beach", "dew line rd", ExactSite))
  
CBAY_2024_siteslook <- CBAY2024_sites %>%
  group_by(ExactSite) %>%
  summarise(countsites = n_distinct(FieldID)) 
  summarise(countsites = n(FieldID)) 
  
ggplot(CBAY_2024_siteslook, aes(y= countsites, x= ExactSite)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) # we have 5 different sites after combing

rm(CBAY_2024_siteslook)

#making a graph for the different black fly species at different sites 

#counting the number of each species at each site

divsitesCBAY20242 <- CBAY2024_sites %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, ExactSite, FieldID, Family, Genus, Species)%>%
  filter(Family == "Simuliidae") %>%
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum*", Species)) %>% #there species were changed based on information taken from NCBI and the phylogetic analysis 
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare*", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium sp. Unclassified", Species)) %>%
  group_by(ExactSite, Species) %>%
  summarise(Speciesnum = n()) %>%
  ungroup()

divsamplelocCBAY2024 <- CBAY2024_sites %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, ExactSite, FieldID, Family, Genus, Species) %>%
  filter(Family == "Simuliidae") %>%
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum", Species)) %>%
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium sp. Unclassified", Species)) %>%
  group_by(ExactSite, FieldID) %>%
  summarise(Speciessum = n_distinct(Species), .groups = "drop")

#plotting this data 

barbfspeciesbysite <-ggplot(divsitesCBAY20242, aes(fill = Species, y= Speciesnum, x= ExactSite)) + 
  geom_bar(position='stack', stat='identity')+
  scale_x_discrete(labels = c("1st Lake"="First Lake", "dew line rd"="Dew Line Road", 
                              "freshwater river"="Freshwater River", "kitigak river"="Kitigak River", 
                              "long pt creek"="Long Point Creek")) +
  scale_fill_manual(values = c("#D62728",  
                               "#F0E442",
                               "#FF5733",  
                               "#CC79A7",
                               "#E69F00",  
                               "#0072B2",  
                               "#009E73",
                               "#FF69B4", "#8A2BE2", "#7AC5CD")) +  
  theme_bw(base_size = 30)+
  xlab("Sampling Location")+
  ylab("Number of Individuals")+
  labs(fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black", margin = margin(t = 5)), axis.title.x = element_text(size = 12),      # X-axis title font size
        axis.title.y = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        legend.title = element_text(size = 12),      # Legend title font size
        legend.text = element_text(size = 10))

barbfspeciesbysite

ggsave("../plots/barsrbf.png", plot = barbfspeciesbysite, width = 9, height = 5, dpi = 300)


#creating a matrix of what black fly species are at which site based on P/A

#turning the data into a P/A matrix

presence_matrix <- divsitesCBAY20242 %>%
  mutate(Speciesnum = ifelse(Speciesnum > 0, 1, 0)) %>%
  distinct() %>%  # remove duplicates
  pivot_wider(names_from = ExactSite, values_from = Speciesnum, values_fill = 0) 

#converting this P/A into the long format

mat_long <- presence_matrix %>%
  as.data.frame() %>%
  pivot_longer(-Species, names_to = "Site", values_to = "Presence") %>%
  mutate(Presence = factor(Presence, levels = c(0,1)))

#ploting the tile plot

tileplot <-ggplot(mat_long, aes(x = Site, y = Species, fill = Presence)) +
  scale_x_discrete(labels = c("1st Lake"="FIL", "dew line rd"="DEW", 
                              "freshwater river"="FWR", "Grenier lake"="GNL", 
                              "long pt creek"="LPC")) +
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

#looking at statistical difference in species richness between sites 

#create data set by coutning the number of each species in each sample

divsamplelocCBAY2024 <- CBAY2024_sites %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, ExactSite, FieldID, Family, Genus, Species) %>%
  filter(Family == "Simuliidae") %>%
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum", Species)) %>%
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium sp. Unclassified", Species)) %>%
  group_by(ExactSite, FieldID) %>%
  summarise(Speciessum = n_distinct(Species), .groups = "drop")

#exploring data for normality and equality of variance

leveneTest(Speciessum ~ ExactSite, data = divsamplelocCBAY2024) #p-value greater than 0.05 - vairances equal 
shapiro.test(divsamplelocCBAY2024$Speciessum) #p-value less than 0.05 - data not normal

#running a non-parmetric test because the data did not meet the above requirements

modeluniqspsample <- kruskal.test(Speciessum ~ ExactSite, data = divsamplelocCBAY2024)

modeluniqspsample



#### PART 3 - Analysis of Sampling Types ----

#looking at different sampling methods to make sure the data makes sense 

CBAY2024_samplingmeth <- CBAY2024_sampledata %>%
  select(SamplingProtocol) %>%
  count(SamplingProtocol, name = "countsamplingtype")

ggplot(CBAY2024_samplingmeth, aes(y= countsamplingtype, x= SamplingProtocol)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))

#creating the data set of the abundence with each sampling type 

divsamplingtypeCBAY2024 <- CBAY2024_sampledata %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, SamplingProtocol, FieldID, Family, Genus, Species)%>%
  filter(Family == "Simuliidae") %>%      
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum*", Species)) %>%
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare*", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium subpusillum*", Species)) %>%
  group_by(SamplingProtocol, Species) %>%
  summarise(Speciessum = n())

view(divsamplingtypeCBAY2024)

#exploring normality and equality of variance 

leveneTest(Speciessum ~ SamplingProtocol, data = divsamplingtypeCBAY2024) #p-value greater than 0.05 - vairances equal 
shapiro.test(divsamplingtypeCBAY2024$Speciessum) #p-value less than 0.05 - data not normal

#running the non-parametric analysis 

modeluniqspsample <- kruskal.test(Speciessum ~ SamplingProtocol, data = divsamplingtypeCBAY2024)

modeluniqspsample

#looking at the average species richness 

summary_divsamplingtypeCBAY2024 <- divsamplingtypeCBAY2024 %>%
  group_by(SamplingProtocol) %>%
  summarise(
    n = n(),
    mean_richness = mean(UniqueSpecies),
    se_richness = sd(UniqueSpecies) / sqrt(n()),
    .groups = "drop")

#exploring this visually 

ggplot(summary_divsamplingtypeCBAY2024, aes(y= mean_richness, x= SamplingProtocol)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))

#comparing the Species richness at each sampling type 

divsamplingtypeCBAY20242 <- CBAY2024_sampledata %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, SamplingProtocol, FieldID, Family, Genus, Species)%>%
  filter(Family == "Simuliidae") %>%
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum*", Species)) %>%
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare*", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium sp. Unclassified", Species)) %>%
  group_by(SamplingProtocol, FieldID) %>%
  summarise(Speciessum = n_distinct(Species), .groups = "drop") %>%
  ungroup()

leveneTest(Speciessum ~ SamplingProtocol, data = divsamplingtypeCBAY20242) #p-value greater than 0.05 - vairances equal 
shapiro.test(divsamplingtypeCBAY20242$Speciessum) #p-value less than 0.05 - data not normal

modeluniqspsample <- kruskal.test(Speciessum ~ SamplingProtocol, data = divsamplingtypeCBAY20242)

modeluniqspsample

#looking at the average for each 

summary_divsamplingtypeCBAY20242 <- divsamplingtypeCBAY20242 %>%
  group_by(SamplingProtocol) %>%
  summarise(
    n = n(),
    mean_richness = mean(Speciessum),
    se_richness = sd(Speciessum) / sqrt(n()),
    .groups = "drop")

#visualizing this 

ggplot(summary_divsamplingtypeCBAY20242, aes(y= mean_richness, x= SamplingProtocol)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))

#### PART 4 - Analysis with river vs terrestrial metadata ----

#combing data with this metadata 

divhabtypeCBAY20242 <- covarraw %>%
  inner_join(divsamplelocCBAY2024, join_by(Sample == FieldID)) 

#testing model assumptions

leveneTest(Speciessum ~ Habitat.type, data = divhabtypeCBAY20242) #p-value greater than 0.05 - vairances equal 
shapiro.test(divhabtypeCBAY20242$Speciessum) #p-value less than 0.05 - data not normal

#running non-parametric test

modeluniqspsample <- kruskal.test(Speciessum ~ Habitat.type, data = divhabtypeCBAY20242)

modeluniqspsample

#looking at summary data 

summary_divhabtypeCBAY20242 <- divhabtypeCBAY20242 %>%
  group_by(Habitat.type) %>%
  summarise(
    n = n(),
    Average = mean(Speciessum),
    SE = sd(Speciessum) / sqrt(n()),
    .groups = "drop") %>%
  mutate(type = "Species Richness") 

#looking at the abundence at each location grouping 

abdunhabtypeCBAY20242 <- covarraw %>%
  inner_join(cbay2024_abundence, join_by(Sample == Sample)) %>%
  mutate(across(c(mosquito_abun, blackfly_abun), ~ ifelse(is.na(.) | . == "", "0", .))) 

#testing model assumptions for balck flies 

leveneTest(blackfly_abun ~ Habitat.type, data = abdunhabtypeCBAY20242) #p-value greater than 0.05 - vairances equal 
shapiro.test(abdunhabtypeCBAY20242$blackfly_abun) #p-value less than 0.05 - data not normal

#running non-parametric test for black flies

modeluniqspsample <- kruskal.test(blackfly_abun ~ Habitat.type, data = abdunhabtypeCBAY20242)

modeluniqspsample

#summary stats for black flies

summary_abdunhabtypeCBAY20242 <- abdunhabtypeCBAY20242 %>%
  group_by(Habitat.type) %>%
  summarise(
    n = n(),
    Average = mean(blackfly_abun),
    SE = sd(blackfly_abun) / sqrt(n()),
    .groups = "drop") %>%
  mutate(type = "Abundence") %>%
  filter(Habitat.type != "Pond")

summarydataforgraph <- bind_rows(summary_divhabtypeCBAY20242, summary_abdunhabtypeCBAY20242)

#graph for abundance and species richness of black flies 

barsumhabdata <- ggplot(summarydataforgraph, aes(y= Average, x= Habitat.type)) +
  geom_col(fill = "#CC79A7") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_errorbar(
    data = summarydataforgraph,
    aes(x = Habitat.type, ymin = Average - SE, ymax = Average + SE),
    width = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x="Habitat Type") +
  theme( panel.border = element_rect(linewidth = 1), 
         axis.text.x = element_text(angle = 0, size = 10, color = "black", margin = margin(t = 5)), 
         axis.title.x = element_text(size = 12, face = "bold"),      # X-axis title font size
        axis.title.y = element_text(size = 12, face = "bold"), axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        legend.title = element_text(size = 12),      # Legend title font size
        legend.text = element_text(size = 10), 
        strip.background = element_rect(color = "black", linewidth = 1, fill = "white"),
        strip.text = element_text(face = "bold", size = "12")) +
  facet_wrap(~type, scales = "free_y")
barsumhabdata
ggsave("../plots/barsumhab.png", plot = barsumhabdata, width = 5, height = 3, dpi = 300)

#making a species composition grpah for the black flies at terrestrial and river locations

divhabtypeCBAY2024 <- CBAY2024_sites %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, ExactSite, FieldID, Family, Genus, Species)%>%
  filter(Family == "Simuliidae") %>%
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum*", Species)) %>%
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare*", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium sp. Unclassified", Species)) %>%
  inner_join(covarraw, join_by(FieldID == Sample)) %>%
  group_by(Habitat.type, Species) %>%
  summarise(Speciesnum = n()) %>%
  ungroup() %>%
  group_by(Habitat.type) %>%
  arrange(Speciesnum, .by_group = TRUE) %>%
  mutate(Species = factor(Species, levels = unique(Species)))

divhabtypeCBAY2024$Species <- factor(divhabtypeCBAY2024$Species , 
                                     levels = rev(c("Metacnephia borealis", "Cnephia eremites", 
                                                "Simulium sp. Unclassified", "Simulium subpusillum","Simulium noelleri","Simulium baffinense",
                                                "Simulium decimatum","Simulium vulgare*","Simulium decorum*","Stegopterna emergens")))


#plotting grpah 

barbfspeciesbyhabtype <-ggplot(divhabtypeCBAY2024, aes(fill = Species, y= Speciesnum, x= Habitat.type)) + 
  geom_bar(position='stack', stat='identity')+
  scale_fill_manual(values = c("#009E73",
                               "#7AC5CD","#0072B2","#FF69B4","#8A2BE2","#E69F00","#CC79A7","#FF5733","#F0E442","#D62728")) +  
  theme_bw(base_size = 30)+
  xlab("Habitat Type")+
  ylab("Number of Individuals")+
  labs(fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black", margin = margin(t = 5)), axis.title.x = element_text(size = 12),      # X-axis title font size
        axis.title.y = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        legend.title = element_text(size = 12),      # Legend title font size
        legend.text = element_text(size = 10))

barbfspeciesbyhabtype

ggsave("../plots/barshabtype.png", plot = barbfspeciesbyhabtype, width = 6, height = 3, dpi = 300)


#### PART 5 - Preparing data for Multidimentional analysis ----

#conduct vegan NMDS analysis 

#create a data frame which works for vegan 

species_bf_ht <-  CBAY2024_sites %>%
  inner_join(CBAYCOI_Domcontigs, join_by(SampleID == Sample)) %>%
  select(SampleID, FieldID, Family, Genus, Species, SamplingProtocol, ExactSite) %>%
  filter(Family == "Simuliidae") %>%
  mutate(Species = ifelse(SampleID == "KBIMP_003_F2", "Simulium decorum", Species)) %>%
  mutate(Species = ifelse(SampleID == "KBIMP-_005_D10", "Simulium vulgare", Species)) %>%
  mutate(Species = ifelse(Species == "unknown", "Simulium sp. Unclassified", Species)) %>%
  count(FieldID, Species, ExactSite, name = "Abundance") %>%
  inner_join(covarraw, join_by(FieldID == Sample))

species_matrix_bfht <- species_bf_ht %>%
  select(FieldID, Species, Abundance) %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "FieldID")

#create the environmental data 

Environmental_bfht <- species_bf_ht %>%
  select(FieldID, Habitat.type, ExactSite) %>%
  distinct(FieldID, Habitat.type, ExactSite)

#run NMDS

NMDS2023 <- metaMDS(species_matrix_bfht, k=3)
en2023 <- envfit(NMDS2023, Environmental_bfht, permutations = 999, na.rm = TRUE)

nmds_scores_2023 <- as.data.frame(vegan::scores(NMDS2023, display = "sites"))
nmds_scores_2023$site <- rownames(nmds_scores_2023)
nmds_scores_2023 <- left_join(nmds_scores_2023, Environmental_bfht, join_by(site == FieldID))

NMDS2023$stress

# Create a function to calculate hulls for each group
calculate_hull <- function(data) {
  data[chull(data$NMDS1, data$NMDS2), ]  # Select convex hull points
}

# Apply the function to each group to create the polygons on the ggplot nmds 
hulls_2023 <- nmds_scores_2023 %>%
  group_by(Habitat.type) %>%
  group_split() %>%
  lapply(calculate_hull) %>%
  bind_rows()


#making the ggplot nmds plot

nmdsplot2023 <- ggplot(data = nmds_scores_2023, aes(x = NMDS1, y = NMDS2)) + 
  geom_polygon(data = hulls_2023, aes(x = NMDS1, y = NMDS2, group = Habitat.type, fill = Habitat.type), 
               alpha = 0.2, colour = "black") +
  geom_point(data = nmds_scores_2023, aes(colour = ExactSite), size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("#7AC5CD","#E69F00")) +  # Map fill colors to sites
  scale_colour_manual(values = c("#00D400","#FF2A7F","#008066","#FF6600"),
                                 labels = c("Grenier lake"="Grenier lake", "dew line rd"="Dew Line Road", 
                                            "freshwater river"="Freshwater River", 
                                            "long pt creek"="Long Point Creek")) +  # Map point colors to ponds
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) + 
  labs(colour = "Sampling location", fill = "Habitat type")
nmdsplot2023
ggsave("../plots/nmdsplothabtypeandlloc.png", nmdsplot2023, width = 4, height = 3, dpi = 300)

anosim(species_matrix_bfht, Environmental_bfht$Habitat.type, permutations = 999, distance = "bray", strata = NULL)
anosim(species_matrix_bfht, Environmental_bfht$ExactSite, permutations = 999, distance = "bray", strata = NULL)

#### - PART 6 - Beta diveristy analysis ----

#looking at the betadiveristy for the different sites

#Creating the data set by combining the vegan data

speciesbetabf <- species_bf_ht %>%
  count(ExactSite, Species, name = "Abundance") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "ExactSite") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

#run beatdiverity analysis 

betacbaybf <- beta.pair(speciesbetabf, index.family = "sorensen")

#explore the beta diveristy results 

dataframebeta <- betacbaybf %>%
  as.data.frame() %>%
  # Convert dist columns to vectors
  transmute(
    beta.sne = as.vector(beta.sne),
    beta.sim = as.vector(beta.sim)
  ) %>%
  # Pivot to long format
  pivot_longer(
    cols = everything(),
    names_to = "type",
    values_to = "Value"
  ) %>%
  mutate(type = case_when(
    type == "beta.sne" ~ "Nestedness.cbay",
    type == "beta.sim" ~ "Turnover.cbay",
    TRUE ~ type  # Default to keeping the original value if no match
  )) 


#create summary data of beta diversity 

summary_betadivcbaybf <- dataframebeta %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean_richness = mean(Value),
    se_richness = sd(Value) / sqrt(n()),
    .groups = "drop")

summary_betadivcbaybf 

#testing the betadiversity values to make sure that they meet the model assumptions

leveneTest(Value ~ type, data = dataframebeta) #p-value greater than 0.05 - vairances equal 
shapiro.test(dataframebeta$Value) #p-value less than 0.05 - data not normal

#running non-parametric T-test

wilcox.test(Value ~ type, data = dataframebeta)

#### PART 7 - analyzing abundence data based on site and between mosquitoes and black flies ----

cbay2024_abundence <- read.csv(file = "../raw-data/cbay2024_abundence.csv")

view(cbay2024_abundence)

#fixing site locations as was done above 

cbay2024_abundence <- cbay2024_abundence %>%
  filter(!is.na(Site) & Site != "") %>%
  mutate(blackfly_abun = ifelse(is.na(blackfly_abun) | blackfly_abun == "", 0, blackfly_abun)) %>%
  mutate(mosquito_abun = ifelse(is.na(mosquito_abun) | mosquito_abun == "", 0, mosquito_abun)) %>%
  mutate(Site = gsub("road", "rd", tolower(Site))) %>%
  mutate(Site = gsub(" resample", "", Site)) %>%
  mutate(Site = gsub("- site 1", "", Site)) %>%
  mutate(Site = gsub("stream ", "", Site)) %>%
  mutate(Site= gsub("west arm/gravel pit near sea", "gravel pit", Site)) %>%
  mutate(Site= gsub("gravel pit/ west arm", "gravel pit", Site)) %>%
  mutate(Site= gsub("west arm", "gravel pit", Site)) %>%
  mutate(Site= gsub("river between", "between", Site)) %>%
  mutate(Site= gsub("grenier lake", "Grenier lake", Site)) %>%
  mutate(Site= gsub("first", "1st", Site)) %>%
  mutate(Site = gsub("\\s+", " ", Site)) %>%
  mutate(Site= gsub("1st creek/", "", Site)) %>%
  mutate(Site= gsub("kitigak river", "Grenier lake", Site)) %>%
  mutate(Site= gsub("2nd culvert pelly rd", "Grenier lake", Site)) %>%
  mutate(Site= gsub("into grenier lake", "Grenier lake", Site)) %>%
  mutate(Site= gsub("1st lake", "Grenier lake", Site)) %>%
  mutate(Site= gsub("between 1st lake and grenier", "Grenier lake", Site)) %>%
  mutate(Site= gsub("inland area", "freshwater river", Site)) %>%
  mutate(Site= gsub("creek south of pelly rd", "freshwater river", Site)) %>%
  mutate(Site= gsub("s of pelly rd near swimming hole", "freshwater river by picnic area", Site)) %>%
  mutate(Site= gsub("swimming hole south of pelly rd", "freshwater river by picnic area", Site)) %>%
  mutate(Site = str_replace(Site, "\\s+$", "")) %>%
  mutate(Site= gsub("freshwater river by picnic area", "freshwater river", Site)) %>%
  mutate(Site= gsub("south of pelly rd", "freshwater river", Site)) %>%
  mutate(Site= gsub("long pt creek near ravine", "long pt creek", Site)) %>%
  mutate(Site= gsub("the ravine", "long pt creek", Site)) %>%
  mutate(Site= gsub("1st creek", "long pt creek", Site)) %>%
  mutate(Site= gsub("past old town", "freshwater river", Site)) %>%
  mutate(Site= gsub("dew line rd freshwater creek", "freshwater river", Site)) %>%
  mutate(Site= gsub("freshwater river pelly rd", "freshwater river", Site)) %>%
  mutate(Site= gsub("freshwater river by swimming hole", "freshwater river", Site)) %>%
  mutate(Site= gsub("2nd culvert pelly rd", "Grenier lake", Site)) %>%
  mutate(Site= gsub("creek by long point beach", "dew line rd", Site)) %>%
  mutate(Site= gsub("pond near 1st lake", "Grenier lake", Site)) %>%
  mutate(Site= gsub("gravel pit/ gravel pit", "gravel pit", Site)) %>%
  mutate(Site= gsub("pelly rd", "freshwater river", Site)) %>%
  mutate(Site= gsub("near dew line station near seasonal creek", "dew line rd", Site)) %>%
  mutate(Site= gsub("augustus hills", "long pt creek", Site)) %>%
  filter(Site != "airport rd pond" & Site != "kuugaq pond") %>%
  filter(Sample != "" & !is.na(Sample)) %>%
  mutate(across(c(mosquito_abun, blackfly_abun), ~ ifelse(is.na(.) | . == "", "0", .))) %>%
  mutate(blackfly_abun = as.numeric(blackfly_abun)) %>%
  mutate(mosquito_abun = as.numeric(mosquito_abun)) 

#pivoting data into the correct format

cbay2024_abundence2 <- cbay2024_abundence %>%
pivot_longer(
    cols = c(blackfly_abun, mosquito_abun),      # columns to pivot
    names_to = "abun_Type",         # new column for names
    values_to = "abundence"              # new column for values
  ) %>%
  mutate(across(c(abundence), ~replace_na(., 0)))
  
#looking at data to make sure it makes sense 

CBAY_2024_siteslookabun <- cbay2024_abundence %>%
  group_by(Site) %>%
  summarise(total_bfabundance = sum(blackfly_abun, na.rm = TRUE))

#visualizing 

ggplot(CBAY_2024_siteslookabun, aes(y= total_bfabundance, x= Site)) +
  geom_col(fill = "skyblue3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))


#analysing difference between black flies and mosquitoes

#testing model assumptions 

leveneTest(abundence ~ abun_Type, data = cbay2024_abundence2) #p-value greater than 0.05 - vairances equal 
shapiro.test(cbay2024_abundence2$abundence) #p-value less than 0.05 - data not normal

#preforming non-parmtric t-test comparing mosquitoes and black flies

wilcox.test(abundence ~ abun_Type, data = cbay2024_abundence2)

#finding summary stats for both 

summary_abungroupCBAY2024 <- cbay2024_abundence2 %>%
  group_by(abun_Type) %>%
  summarise(
    n = n(),
    mean_abun = mean(abundence),
    se_abun = sd(abundence) / sqrt(n()),
    .groups = "drop") 

#Analysing black flies based on site

#Visualizing normality 
boxplot(abun_type ~ Site, data = cbay2024_abundence2,
        col = "red",
        border = "black", 
        las = 2) #make x-axis labels perpendicular

#testing model assumptions
leveneTest(blackfly_abun ~ Site, data = cbay2024_abundence) #p-value greater than 0.05 - vairances equal 
shapiro.test(cbay2024_abundence$blackfly_abun) #p-value less than 0.05 - data not normal

#running non-parametric model 
modeluniqspsample <- kruskal.test(blackfly_abun ~ Site, data = cbay2024_abundence)

modeluniqspsample

#testing model assumptions for mosquitoes
leveneTest(mosquito_abun ~ Site, data = cbay2024_abundence) #p-value greater than 0.05 - vairances equal 
shapiro.test(cbay2024_abundence$mosquito_abun) #p-value less than 0.05 - data not normal

#running non-parametric model
modeluniqspsample <- kruskal.test(mosquito_abun ~ Site, data = cbay2024_abundence)

modeluniqspsample

#finding summary stats based on site
summary_abunsitesCBAY20242 <- cbay2024_abundence %>%
  group_by(Site) %>%
  summarise(
    n = n(),
    mean_bfabun = mean(blackfly_abun),
    mean_mosabun = mean(mosquito_abun),
    se_bfabun = sd(blackfly_abun) / sqrt(n()),
    se_mosabun = sd(mosquito_abun) / sqrt(n()),
    .groups = "drop") %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "Type",
               values_to = "mean") %>%
  pivot_longer(cols = starts_with("se"),
               names_to = "Typese",
               values_to = "se") %>%
  mutate(Type = gsub("mean_", "", Type))


