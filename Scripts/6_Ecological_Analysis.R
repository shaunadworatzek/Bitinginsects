
#To run this data set you need the meatdata from 2024, the fixed site names from 
#2024, the metadata from both places in 2025, the updates species data sets from 
#both years (KBIMP2025_updatedspecies and KBIMP2025_updatedspecies), the data set
#from schafer from 2012 (sr_2012), and the vector dataset from this time 
#this will output multiple plots including the rarefaction curves, vector change plot
#total species richness, and multidimensonal sclaing and well as mulitple statistical tests. 

#### PART 1 - Loading packages, Openeing files, Running functions ----

##### Packages #####


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
library(maps)
library(ggmap)
library(sf)
library(rgbif)
library(rnaturalearth)
library(paletteer)
library(gt)
library(cancensus)


##### Files #####

CBAY2025_metadata <- read_csv(file = "raw-data2/CBAY2025_metadata.csv")
KGLTK2025_metadata <- read_csv(file = "raw-data2/KGLTK2025_metadata.csv")
KBIMP_combined <- read_tsv(file = "processed-data/KBIMP_updatedspecies.tsv")
vector_change <- read_csv(file = "raw-data2/vector_change.csv")
kbimp2024_sampledata_clean <- read_csv(file = "processed-data/kbimp2024_sampledata_clean.csv")
nonbiting_species <- read.csv(file = "raw-data2/non_biting_species.csv")
expectedbf <- read.csv(file = "raw-data2/expectedbf.csv")
NBP <- read.csv(file = "raw-data2/NorthernBiodiversityProgram - Simuliidae P.Schaefer.csv")
NBP2 <- read.csv(file = "raw-data2/ROM_NBP_Simuliidae_2.csv")
NBP3 <- read.csv(file = "raw-data2/ROM-NPB Simuliidae3.csv")
KBIMP2024 <- read_tsv(file = "processed-data/KBIMP2024_clean.tsv")

##### Packages #####

# Create a function to calculate hulls for each group
calculate_hull <- function(data) {
  data[chull(data$NMDS1, data$NMDS2), ]  # Select convex hull points
}

#function to fix coordinates

clean_coords <- function(x) {
  x %>%
    str_replace_all("\\.\\s+", ".") %>%  # fix ". 4827" → ".4827"
    str_replace("^([0-9]+\\s+[0-9]+)\\s+([0-9]+)$", "\\1.\\2")  # add missing period
}



#### PART 2 - Preparing the site metadata for analysis ----

##### preparing the 2025 data #####

#creating month field and selecting required fields in CBAY 

CBAY2025_metadata_clean <- CBAY2025_metadata %>%
  mutate(Date_parsed = 
          parse_date_time(`Date collected`,
          orders = c("dmy")),
         Month = month(Date_parsed, label = TRUE)) %>%
  select(Month, Sample, Site,  `Sample type collection method`,
         Lon, Lat, `Habitat type`,  `Mosquito Head Abundance`, 
         `Blackfly Head Abundance`) %>%
  filter(!is.na(`Mosquito Head Abundance`)) 

#creating month field, reformting lat/lon, and selecting required fields in Kugluktuk 

KGLTK2025_metadata_clean <- KGLTK2025_metadata %>%
  mutate(Date_parsed = 
           parse_date_time(`Date set`,
            orders = c("ymd", "mdy", "dmy", "dm")),
         Month = month(Date_parsed, label = TRUE)) %>% #some of the dates are parsed incorrectly but we doubled checked that the correct months were pulled. 
  filter(!is.na(`Mosquito Head Abundance`)) %>%
  select(Month, Sample, Site, `Sample type collection method`,
         Lon, Lat, `Habitat type`,  `Mosquito Head Abundance`, 
         `Blackfly Head Abundance`) %>%
  separate(Lat, into = c("Latdeg", "Latmin"), sep = " ", 
           convert = TRUE) %>%
  mutate(Lat = Latdeg + (Latmin / 60)) %>%
  separate(Lon, into = c("Londeg", "Lonmin"), sep = " ", 
           convert = TRUE) %>%
  mutate(Lon = -(Londeg + (Lonmin / 60))) %>%
  select(-Latmin, -Lonmin, -Latdeg, -Londeg)

#combining KGLTKand CBAY data and creating Sector field

KBIMP2025_metadata <- CBAY2025_metadata_clean %>% 
  rbind(KGLTK2025_metadata_clean) %>%
  mutate(Sector = str_extract(Sample, "^[A-Za-z]+")) %>%
  mutate(Year = 2025) %>%
  dplyr::rename(SamplingMethod = `Sample type collection method`) 

##### preparing the 2024 data #####

KBIMP2024_metadata <- kbimp2024_sampledata_clean %>%
  dplyr::rename(SamplingMethod = SamplingProtocol) %>%
  dplyr::rename(Sample = FieldID) %>%
  select(Sample, SamplingMethod, Sector, Month, Lon, Lat) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  separate(Lat, into = c("Latdeg", "Latmin"), sep = " ",
           convert = TRUE, remove = FALSE) %>%
  mutate(Lat = case_when(Sector == "KGLTK" ~ 
                           Latdeg + (Latmin / 60), 
    TRUE ~ as.numeric(Lat))) %>%
  separate(Lon, into = c("Londeg", "Lonmin"), sep = " ", 
           convert = TRUE, remove = FALSE) %>%
  mutate(Lon = case_when(Sector == "KGLTK" ~ 
           -(Londeg + (Lonmin / 60)), TRUE ~ as.numeric(Lon))) %>%
  select(-Latmin, -Lonmin, -Latdeg, -Londeg)

##### combining 2024 and 2025 data into one metadata file #####

meta_data <- KBIMP2024_metadata %>%
  mutate(Year = 2024) %>%
  bind_rows(KBIMP2025_metadata) %>%
  mutate(SamplingMethod = case_when(
    grepl("malaise", SamplingMethod, ignore.case = TRUE) ~ "Malaise Trap",
    grepl("sweep", SamplingMethod, ignore.case = TRUE) ~ "Sweep Net",
    grepl("aspirator|people", SamplingMethod, ignore.case = TRUE) ~ "Aspirator",
    TRUE ~ SamplingMethod)) %>%
  filter(SamplingMethod %in% 
           c("Malaise Trap", "Sweep Net", "Aspirator")) %>%
  select(Sample, SamplingMethod, Sector, Month, Lon, Lat, Year) 

#removing the outgroup from the data
KBIMP_combined <- KBIMP_combined %>%
  mutate(Sample = str_remove(Sample, "_.*")) %>%
  filter(!is.na(Species)) %>% 
  distinct(across(-update_flag), .keep_all = TRUE) %>%
  filter(Sample != "Outgroup")

#### PART 3 - Creating Sample/Site number and PA tables  ----

##### 2024 #####

#number of samples 

tabledata2 <- KBIMP2024_metadata %>%
  left_join(KBIMP_combined, join_by(Sample == Sample)) %>%
  select(Sample, Sector, Family) %>%
  group_by(Sample, Sector) %>%
  summarise(has_mosquito = any(Family == "Culicidae", na.rm = TRUE),
    has_blackfly = any(Family == "Simuliidae", na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    has_mosquito = as.integer(has_mosquito),
    has_blackfly = as.integer(has_blackfly)) %>%
 mutate(category = case_when(
    has_mosquito != 0 & has_blackfly != 0  ~ "Both",
    has_mosquito != 0 & has_blackfly == 0 ~ "Mosquito only",
    has_mosquito == 0 & has_blackfly != 0 ~ "Black fly only",
    has_mosquito == 0 & has_blackfly == 0 ~ "Neither")) %>%
  group_by(Sector, category) %>%
  summarise(value = n_distinct(Sample)) %>%
  ungroup() 

#total number of samples

total_row <- tabledata2 %>%
  group_by(Sector) %>%
    summarise(
      category = "Total",
      value = sum(value)) %>%
    ungroup()

#combing into one for 2024 data
  
sample_table_2024 <- bind_rows(tabledata2, total_row) %>%
  mutate(year = 2024)
  
#number of individuals 

individual_table_2024 <- KBIMP2024 %>%
  mutate(Sector = 
           str_extract(Sample, "^[A-Za-z]+")) %>%
  select(Sample, Sector, Family) %>%
  group_by(Sector) %>%
  summarise(Blackflies = sum(Family == "Simuliidae", 
                         na.rm = TRUE),
            Mosquitoes = sum(Family == "Culicidae", 
                         na.rm = TRUE)) %>%
  pivot_longer(cols = c("Blackflies", "Mosquitoes"), values_to = "value", names_to = "category") %>%
  mutate(year = 2024)

##### 2025 #####

# number of individuals 

individuals_table2025 <- KBIMP2025_metadata %>%
  mutate(Sector = str_extract(Sample, "^[A-Za-z]+")) %>%
  select(Sector, `Mosquito Head Abundance`,`Blackfly Head Abundance`) %>%
  dplyr::rename("Mosquitoes" = "Mosquito Head Abundance", 
         "Blackflies" = "Blackfly Head Abundance") %>%
  pivot_longer(cols = c("Mosquitoes", "Blackflies"), values_to = "value", names_to = "category") %>%
  group_by(Sector, category) %>%
summarise(
  value = sum(value)) %>%
  ungroup() %>%
  mutate(year = 2025)

#number of samples 
  
samples_table2025 <- KBIMP2025_metadata %>%
  mutate(Sector = str_extract(Sample, "^[A-Za-z]+")) %>%
  left_join(KBIMP_combined, join_by(Sample == Sample)) %>%
  select(Sample, Sector, Family) %>%
  group_by(Sample, Sector) %>%
  summarise(has_mosquito = any(Family == "Culicidae", na.rm = TRUE),
            has_blackfly = any(Family == "Simuliidae", na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    has_mosquito = as.integer(has_mosquito),
    has_blackfly = as.integer(has_blackfly)) %>%
  mutate(category = case_when(
    has_mosquito != 0 & has_blackfly != 0  ~ "Both",
    has_mosquito != 0 & has_blackfly == 0 ~ "Mosquito only",
    has_mosquito == 0 & has_blackfly != 0 ~ "Black fly only",
    has_mosquito == 0 & has_blackfly == 0 ~ "Neither")) %>%
  group_by(Sector, category) %>%
  summarise(value = n_distinct(Sample)) %>%
  ungroup() 

#total number of samples 

total_row2025 <- samples_table2025 %>%
  group_by(Sector) %>%
  summarise(
    category = "Total",
    value = sum(value)) %>%
  ungroup()

#combining into one for 2025

samples_table2025 <- bind_rows(samples_table2025, total_row2025) %>%
  mutate(year = 2025) 

##### individuals table final for both years #####

individuals_table <- 
  bind_rows(individual_table_2024, individuals_table2025) %>%
  pivot_wider(names_from = "Sector", values_from = "value") %>%
  dplyr::rename("Cambridge Bay\n(Iqaluktuuttiaq)\n" = "CBAY", 
         "Kugluktuk\n(Qurluqtuq)\n" = "KGLTK")

individuals_gtable <- individuals_table %>%
  gt(rowname_col = "category", groupname_col = "year") %>%
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
    locations = list(cells_column_labels(),
      cells_row_groups())) %>%
  tab_style(style = cell_borders(sides = c("bottom"),
            color = "black", weight = px(2)),
            locations = cells_title()) %>%
    tab_style(style = cell_borders(sides = c("top"),
             color = "black", weight = px(2)),
              locations = cells_footnotes()) %>%
  tab_style(style = cell_borders(
    sides = c("top", "bottom"),
            color = "black", weight = px(2)),
            locations = cells_row_groups()) %>%
  tab_options(data_row.padding = px(5)) %>%
  cols_width("Cambridge Bay\n(Iqaluktuuttiaq)\n" 
             ~ px(300), 
             "Kugluktuk\n(Qurluqtuq)\n" ~ px(300))

gtsave(data = individuals_gtable, 
       filename = "plots/individualstable2.png")

##### samples table final for both years #####

samples_table <- bind_rows(sample_table_2024, samples_table2025) %>%
  pivot_wider(names_from = "Sector", values_from = "value") %>%
  dplyr::rename("Cambridge Bay\n(Iqaluktuuttiaq)\n" = "CBAY", 
         "Kugluktuk\n(Qurluqtuq)\n" = "KGLTK") %>%
  mutate(category = factor(category, 
        levels = c("Total", "Black fly only", 
                   "Mosquito only", "Both", "Neither"))) %>%
  arrange(category)


samples_gtable <- samples_table %>%
  gt(rowname_col = "category", groupname_col = "year") %>%
  tab_header(title = md("")) %>%
  tab_footnote(footnote = md("")) %>%
  cols_align(align = "right",
    columns = everything()) %>%
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
            locations = list(cells_column_labels(),
                             cells_row_groups())) %>%
  tab_style(style = cell_borders(sides = c("bottom"),
            color = "black", weight = px(2)),
            locations = cells_title()) %>%
  tab_style(style = cell_borders(sides = c("top"),
             color = "black", weight = px(2)),
            locations = cells_footnotes()) %>%
  tab_style(style = cell_borders(sides = c("top", "bottom"),
            color = "black", weight = px(2)),
            locations = cells_row_groups()) %>%
  tab_options(data_row.padding = px(5)) %>%
  cols_width("Cambridge Bay\n(Iqaluktuuttiaq)\n" 
             ~ px(300), 
             "Kugluktuk\n(Qurluqtuq)\n" ~ px(300))

gtsave(data = samples_gtable, 
       filename = "plots/samples_gtable.png")


##### PA table for thesis ##### 

PA_gtable <- KBIMP_combined %>%
  mutate(Species = if_else(Species == "Aedes nigripes/impiger", "Aedes nigripes/Aedes impiger", Species)) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  mutate(presence = 1) %>%
  select(-Sample, -Family, -update_flag, -Genus, -Sequence) %>%
  distinct() %>%  #I think I fixed the issue with the code you pointed out here - I get the same result but I think this is the better way of doing this 
  pivot_wider(names_from = "region", 
              values_from = "presence", 
              values_fill = 0) %>%
  dplyr::rename("Cambridge Bay (Iqaluktuuttiaq)" = "CBAY", 
         "Kugluktuk (Qurluqtuq)" = "KGLTK") %>%
  arrange(Species) %>%
  mutate(Family = if_else(Species %in% 
                            c("Aedes communis",
                              "Aedes excrucians",
                              "Aedes nigripes/Aedes impiger",
                              "Aedes punctor/Aedes hexodontus",
                              "Culiseta alaskaensis",
                              "Culiseta inornata", 
                              "Aedes koreicus",
                              "Aedes cataphylla"),
                          "Mosquitoes","Black flies")) %>%
  gt(rowname_col = "Species", groupname_col = "Family") %>%
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
    locations = cells_stub())


gtsave(data = PA_gtable, 
       filename = "plots/PA_gtable.png")


####  PART 4 - Running iNEXT analysis and generating graphs ----

#getting at the number of samples  for iNEXT

SamplesCBAY <- sum(meta_data$Sector == "CBAY",
                   na.rm = TRUE)

SamplesKGLTK <- sum(meta_data$Sector == "KGLTK",
                    na.rm = TRUE)

SamplesJulCBAY <- sum(meta_data$Month == "Jul" &
                        meta_data$Sector == "CBAY",
                      na.rm = TRUE)

SamplesJulKGLTK <- sum(meta_data$Month == "Jul" &
                         meta_data$Sector == "KGLTK",
                       na.rm = TRUE)

SamplesAugCBAY <- sum(meta_data$Month == "Aug" &
                        meta_data$Sector == "CBAY",
                      na.rm = TRUE)

SamplesAugKGLTK <- sum(meta_data$Month == "Aug" &
                         meta_data$Sector == "KGLTK",
                       na.rm = TRUE)

##### iNEXT both places, both families #####  

#reorganizing data into format needed for iNEXT (incidence of each species by region with total number of samples as the first row)

iNEXT <- KBIMP_combined %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Presence") %>% 
  #assigning 1 to all species/ sample combos for presence in a sample
  dplyr::count(region, Species, name = "Incidence") %>%
  #summing all the presence values by region and species to get incidence of each species 
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  #here we are replacing NAs with 0 because some of the species are present in CBAY but not KGLTK and vise versa so this makes the incidence of those species in those regions 0. 
  add_row(Species = "sampling_extent",
          CBAY = SamplesCBAY,
          KGLTK = SamplesKGLTK, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext <- iNEXT(iNEXT, q=0, datatype="incidence_freq")

em.inext$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_total <- as.data.frame(em.inext$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(type = "Total")

em.inext$DataInfo

#creating figure

iNext_total <- ggiNEXT(em.inext)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18)) +
  annotate("text", x = 170, y = 10, label = "99.7%", color = "#000099", size = 5, hjust = 0) +
  annotate("text", x = 65, y = 30, label = "91.8%", color = "#FFC000", size = 5, hjust = 0) +
  ggtitle("Black flies and Mosquitoes Combined")

##### iNEXT both places, just black flies both months #####  

iNEXT_sim <- KBIMP_combined %>%
  filter(Family == "Simuliidae" )%>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Presence") %>% 
  #assigning 1 to all species/ sample combos for presence in a sample
  dplyr::count(region, Species, name = "Incidence") %>%
  #summing all the presence values by region and species to get incidence of each species 
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  add_row(Species = "sampling_extent",
          CBAY = 209,
          KGLTK = 114, .before = 1) %>%
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

#creating figure

iNext_bf <- ggiNEXT(em.inext.sim)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18)) +
  annotate("text", x = 200, y = 7, label = "99.3%", color = "#000099", size = 5) +
  annotate("text", x = 100, y = 23, label = "87.5%", color = "#FFC000", size = 5) +
  ggtitle("Black flies")

ggsave("plots/inext2025bf.png", iNext_bf , width = 6, height = 4, dpi = 300, bg = "transparent")

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
          CBAY = 209,
          KGLTK = 114, .before = 1) %>%
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

#creating figure

iNext_mos <- ggiNEXT(em.inext.cul)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18)) +
  annotate("text", x = 200, y = 1.5, label = "100%", color = "#000099", size = 5) +
  annotate("text", x = 140, y = 8, label = "98.0%", color = "#FFC000", size = 5) +
  ggtitle("Mosquitoes")


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
          CBAY = 120,
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

em.inext.cul.jul$DataInfo

#creating figure

iNext_cul_jul <- ggiNEXT(em.inext.cul.jul)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18)) +
  annotate("text", x = 105, y = 4.5, label = "100%", color = "#000099", size = 5, hjust = 0) +
  annotate("text", x = 50, y = 8, label = "94.3%", color = "#FFC000", size = 5, hjust = 0) +
  annotate("text", x = 200, y = 8, label = "July", angle = 0, size = 6, fontface = "bold")


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
          CBAY = 120,
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

#creating figure

iNext_bf_jul <- ggiNEXT(em.inext.bf.jul)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18)) +
  annotate("text", x = 110, y = 6, label = "95.6%", color = "#000099", size = 5, hjust = 0) +
  annotate("text", x = 30, y = 15, label = "90.3%", color = "#FFC000", size = 5, hjust = 0) +
  annotate("text", x = 200, y = 16, label = "July", angle = 0, size = 6, fontface = "bold")

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
          KGLTK = 56, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext.cul.aug <- iNEXT(iNEXT_cul_aug, q=0, datatype="incidence_freq")

em.inext.cul.aug$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_mos_aug <- as.data.frame(em.inext.cul.aug$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(Family = "Culicidae") %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(Month = "August")

em.inext.cul.aug$DataInfo

#creating figure

iNext_cul_aug <- ggiNEXT(em.inext.cul.aug)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18)) +
  annotate("text", x = 80, y = 2, label = "100%", color = "#000099", size = 5, hjust = 0) +
  annotate("text", x = 30, y = 4, label = "75.9%", color = "#FFC000", size = 5, hjust = 0) +
  annotate("text", x = 135, y = 4.7, label = "August", angle = 0, size = 6, fontface = "bold")


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
          KGLTK = 56, .before = 1) %>%
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

em.inext.bf.aug$DataInfo

#creating figure

iNEXT_bf_aug <- ggiNEXT(em.inext.bf.aug)+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  scale_fill_manual(values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                    breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                    name = "Region") +
  scale_shape_manual( values = c(15, 16, 17, 18)) +
  annotate("text", x = 85, y = 7.5, label = "94.8%", color = "#000099", size = 5, hjust = 0) +
  annotate("text", x = 30, y = 18, label = "80.6%", color = "#FFC000", size = 5, hjust = 0) +
  annotate("text", x = 130, y = 28, label = "August", angle = 0, size = 6, fontface = "bold")


##### Combining plots #####

layout_design <- "
  11
  23
  45
  67"
  

combined_iNextplot <- (iNext_total + iNext_mos + iNext_bf + 
                         iNext_cul_jul + iNext_bf_jul +
                         iNext_cul_aug + iNEXT_bf_aug) +
  plot_layout(design = layout_design, guides = "collect") 

png("plots/combinediNext.png", width = 3500, height = 4000, res = 300)

print(combined_iNextplot)

dev.off()

#### PART 5 - Investigating differences since 2012 ----

##### getting the Scheader data and combing different sources #####

#getting data from gbif



#getting usage key for mosquitoes 

#xC <- name_backbone("Culicidae")
#xC$usageKey

#downloading data key related file

#(gbif_download_mos <- occ_download(
  #pred("country", "CA"),   #selecting for data measured in Canada 
  #pred("taxonKey", xC$usageKey), #for Insecta which I previously set as my usage key 
  #pred("hasGeospatialIssue", FALSE), #make sure none of the data I download has a geospatial issue 
  #pred("hasCoordinate", TRUE),  #make sure all the data has coordinates 
  #pred("occurrenceStatus", "PRESENT"), #make sure that it is occurrence data focusing on presence rather than absence 
  #pred_gte("year", 1900),  #from the year 2000 till the present 
  #pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN"))), #exclude fossil records to try and make sure we are getting live specimens
  #user = user, pwd = pwd, email = email,
  #format = "SIMPLE_CSV")) #log in information and download format

#occ_download_wait(gbif_download_mos)

#opening data into R

my_download_key <- "0037903-260623161305970" 

records_gbif_mos <- occ_download_get(key = my_download_key, path = "plots/") |> 
  occ_download_import()

#filtering data between lat and lon of our sampling in CBAY

records_gbif_CBAY <- records_gbif_mos %>%
  filter(decimalLatitude >= 69  & decimalLatitude <= 69.23) %>%
  filter(decimalLongitude >= -105.8 & decimalLongitude <= -104.77)

#filtering data between lat and lon of our sampling in KGLTK

records_gbif_Kugluktuk <- records_gbif_mos %>%
  filter(decimalLatitude >= 67.79 & decimalLatitude <= 68.15) %>%
  filter(decimalLongitude >= -115.3 & decimalLongitude <= -115)

#getting usage key for black flies

#x <- name_backbone("Simuliidae")
#x$usageKey

#downloading data key related file

#(gbif_download_bf <- occ_download(
 # pred("country", "CA"),   #selecting for data measured in Canada 
  #pred("taxonKey", x$usageKey), #for Insecta which I previously set as my usage key 
#  pred("hasGeospatialIssue", FALSE), #make sure none of the data I download has a geospatial issue 
 # pred("hasCoordinate", TRUE),  #make sure all the data has coordinates 
  #pred("occurrenceStatus", "PRESENT"), #make sure that it is occurrence data focusing on presence rather than absence 
  #pred_gte("year", 1900),  #from the year 2000 till the present 
  #pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN"))), #exclude fossil records to try and make sure we are getting live specimens
  #user = user, pwd = pwd, email = email,
  #format = "SIMPLE_CSV")) #log in information and download format

#occ_download_wait(gbif_download_bf)

#opening data into R

my_download_key <- "0037860-260623161305970" 

records_gbif_bf <- occ_download_get(key = my_download_key, path = "plots/") |> 
  occ_download_import()

#filtering data between lat and lon of our sampling in CBAY

records_gbif_CBAY_bf <- records_gbif_bf %>%
  filter(decimalLatitude >= 69  & decimalLatitude <= 69.23) %>%
  filter(decimalLongitude >= -105.8 & decimalLongitude <= -104.77)

#filtering data between lat and lon of our sampling in KGLTK

records_gbif_Kugluktuk <- records_gbif_mos %>%
  filter(decimalLatitude >= 67.79 & decimalLatitude <= 68.15) %>%
  filter(decimalLongitude >= -115.3 & decimalLongitude <= -115)

#filtering 3 files from the digitization at the ROM

NBP_filtered <- NBP %>%
  filter(Sector %in% c("Cambridge Bay", "Kugluktuk")) %>%
  mutate(Sector = case_when(
    Sector == "Cambridge Bay" ~ "CBAY", 
    Sector == "Kugluktuk"     ~ "KGLTK")) %>%
  mutate(latitude = round(Lat, 2)) %>%
  mutate(longitude = round(Lon, 2)) %>%
  dplyr::rename(family = Family) %>%
  dplyr::rename(species = Species) %>%
  select(family, species, latitude, longitude, Sector) 

NBP2_filtered <- NBP2 %>%
  filter(
    (between(latitude, 69, 69.23) & 
       between(longitude, -105.8, -104.77)) |
      (between(latitude, 67.79, 68.15) & 
         between(longitude, -115.3, -115))) %>% 
  unite(col = "species", c(genus, specificEpithet), sep = " ") %>%
  select(family, species, latitude, longitude) %>%
  mutate(Sector = case_when(
    between(latitude, 69, 69.23) & 
      between(longitude, -105.8, -104.77) ~ "CBAY",
    between(latitude, 67.79, 68.15) 
    & between(longitude, -115.3, -115) ~ "KGLTK")) %>%
  mutate(latitude = round(latitude, 2)) %>%
  mutate(longitude = round(longitude, 2))

NBP3_filtered <- NBP3 %>%
  filter(
    (between(latitude, 69, 69.23) & 
       between(longitude, -105.8, -104.77)) |
      (between(latitude, 67.79, 68.15) & 
         between(longitude, -115.3, -115))) %>% 
  unite(col = "species", c(genus, specificEpithet), sep = " ") %>%
  select(family, species, latitude, longitude) %>%
  mutate(Sector = case_when(
    between(latitude, 69, 69.23) & 
      between(longitude, -105.8, -104.77) ~ "CBAY",
    between(latitude, 67.79, 68.15) 
    & between(longitude, -115.3, -115) ~ "KGLTK")) %>%
  mutate(latitude = round(latitude, 2)) %>%
  mutate(longitude = round(longitude, 2))


#filtering for just the Schaefer data in the gbif bf file 

records_gbif_bf_schaefer <- records_gbif_bf %>%
  filter(
    (between(decimalLatitude, 69, 69.23) & 
       between(decimalLongitude, -105.8, -104.77)) |
      (between(decimalLatitude, 67.79, 68.15) & 
         between(decimalLongitude, -115.3, -115))) %>%
  filter(recordedBy == "P. Schaefer") %>% 
  mutate(Sector = case_when(
    between(decimalLatitude, 69, 69.23) & 
      between(decimalLongitude, -105.8, -104.77) ~ "CBAY",
    between(decimalLatitude, 67.79, 68.15) 
    & between(decimalLongitude, -115.3, -115) ~ "KGLTK")) %>%
  dplyr::rename(latitude = decimalLatitude) %>%
  dplyr::rename(longitude = decimalLongitude) %>%
  mutate(latitude = round(latitude, 2)) %>%
  mutate(longitude = round(longitude, 2)) %>%
  select(family, species, Sector, latitude, longitude)

#filtering for just the Schaefer data in the gbif mosquito file 

records_gbif_mos_schaefer <- records_gbif_mos %>%
  filter(
    (between(decimalLatitude, 69, 69.23) & 
       between(decimalLongitude, -105.8, -104.77)) |
      (between(decimalLatitude, 67.79, 68.15) & 
         between(decimalLongitude, -115.3, -115))) %>%
  filter(recordedBy == "P. Schaefer") %>% 
  mutate(Sector = case_when(
    between(decimalLatitude, 69, 69.23) & 
      between(decimalLongitude, -105.8, -104.77) ~ "CBAY",
    between(decimalLatitude, 67.79, 68.15) 
    & between(decimalLongitude, -115.3, -115) ~ "KGLTK")) %>%
  dplyr::rename(latitude = decimalLatitude) %>%
  dplyr::rename(longitude = decimalLongitude) %>%
  mutate(latitude = round(latitude, 2)) %>%
  mutate(longitude = round(longitude, 2)) %>%
  select(family, species, Sector, latitude, longitude)

#combining the five data files for the Schaefer data 

schaeferdata <- bind_rows(NBP_filtered, NBP2_filtered, NBP3_filtered, 
                          records_gbif_bf_schaefer, 
                          records_gbif_mos_schaefer)

#converting the Schaefer data into a PA format 

schaeferdata_PA <- schaeferdata %>%
  mutate(presence = 1) %>%
  mutate(species = case_when(
    species == "Simulium tuberosum complex" ~ "Simulium tuberosum",
    species == "Metacnephia boreale" ~ "Metacnephia borealis",
    species == "Simulium venustum" ~ "Simulium venustum complex",
    TRUE ~ species)) %>%
  select(-family, -latitude, -longitude) %>%
  distinct() %>% 
  pivot_wider(names_from = "Sector", 
              values_from = "presence", 
              values_fill = 0) %>%
  mutate(species = if_else(species == "Aedes nigripes", "Aedes nigripes/Aedes impiger", species)) %>%
  mutate(species = str_replace(species,  "Aedes hexodontus",
                               "Aedes punctor/Aedes hexodontus")) %>%
  filter(!species == "Aedes impiger") %>%
  filter(!species == "") %>%
  dplyr::rename("Cambridge Bay 2010-2012" = "CBAY", 
         "Kugluktuk 2010-2012" = "KGLTK") 

#producing a finalized table of the lat/lon coordinators rounded 
#to 2 decimals for the Schaefer data 

latlonschaefer_table <- schaeferdata %>%
  select(Sector, latitude, longitude) %>%
  distinct() %>%
  arrange(latitude) %>%
  group_by(Sector) %>%
  mutate(Site_ID = row_number()) %>%
  ungroup() %>%
  gt(rowname_col = "Site_ID", groupname_col = "Sector") %>%
  cols_width(
    stub() ~ px(100),              
    everything() ~ px(200)) %>%
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
  tab_options(data_row.padding = px(5)) 

gtsave(data = latlonschaefer_table, 
       filename = "plots/latlonschaefer.png")


##### venn diagram/ determining which species are different #####

#combining Schaefer data with my data to compare species presence 

venndiagramspecieslist <- KBIMP_combined  %>%
  mutate(Species = if_else(Species == "Aedes nigripes/impiger", 
                           "Aedes nigripes/Aedes impiger", Species)) %>%
  mutate(Species = if_else(Species == "Aedes hexodontus", 
                           "Aedes punctor/Aedes hexodontus", Species)) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  mutate(presence = 1) %>%
  select(-Sample, -Family, -update_flag, -Genus, -Sequence) %>%
  distinct() %>% #I think I fixed the issue with the code you pointed out here - I get the same result but I think this is the better way of doing this 
  pivot_wider(names_from = "region", 
              values_from = "presence", 
              values_fill = 0) %>%
  dplyr::rename("Cambridge Bay 2024-2025" = "CBAY", 
         "Kugluktuk 2024-2025" = "KGLTK") %>%
  full_join(schaeferdata_PA, join_by("Species" == "species")) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  filter(!if_all(-Species, ~ . == 0)) %>%
  mutate(Family = if_else(Species %in% 
                  c("Aedes communis",
                    "Aedes excrucians",
                    "Aedes nigripes/Aedes impiger",
                    "Aedes punctor/Aedes hexodontus",
                    "Culiseta alaskaensis",
                    "Culiseta inornata", 
                    "Aedes koreicus",
                    "Aedes cataphylla"),
        "Mosquitoes","Black flies"))

#just looking at this data in CBAY

venncbay <- venndiagramspecieslist %>%
  select(`Cambridge Bay 2010-2012`, `Cambridge Bay 2024-2025`, Species) %>%
  pivot_longer(-Species, names_to = "sector", values_to = "presence") %>%
  filter(presence == 1) %>%
  select(sector, Species) %>%
  group_by(sector) %>%
  summarise(species_list = list(unique(Species))) %>%
  deframe()

#species in our survey not in the 2010-2012 survey 

setdiff(venncbay$`Cambridge Bay 2024-2025`, 
        venncbay$`Cambridge Bay 2010-2012`)



setdiff(venncbay$`Cambridge Bay 2010-2012`, 
        venncbay$`Cambridge Bay 2024-2025`)

#species in both surveys 

intersect(venncbay$`Cambridge Bay 2010-2012`, 
        venncbay$`Cambridge Bay 2024-2025`)

#just looking at this data in CBAY

vennkug <- venndiagramspecieslist %>%
  select(`Kugluktuk 2024-2025`, 
         `Kugluktuk 2010-2012`, Species) %>%
  pivot_longer(-Species, names_to = "sector", values_to = "presence") %>%
  filter(presence == 1) %>%
  select(sector, Species) %>%
  group_by(sector) %>%
  summarise(species_list = list(unique(Species))) %>%
  deframe()

#species in the 2010-2012 survey which weren't in our survey 

setdiff(vennkug$`Kugluktuk 2010-2012`, 
        vennkug$`Kugluktuk 2024-2025`)

#species in our survey not in the 2010-2012 survey 

setdiff(vennkug$`Kugluktuk 2024-2025`, 
        vennkug$`Kugluktuk 2010-2012`)

#species in both surveys 

intersect(vennkug$`Kugluktuk 2024-2025`, 
        vennkug$`Kugluktuk 2010-2012`)

#### PART 6 - Alpha diversity analysis ----

#create data set by counting the number of each species in each sample

alpha_en_data <- meta_data %>%
  filter(SamplingMethod %in% c("Malaise Trap", "Sweep Net")) 

families <- c("Culicidae", "Simuliidae") 

speciesrich_CBAYvs_KGLTK <- KBIMP_combined %>%
  inner_join(alpha_en_data) %>% #I think I fixed the problem here so we do not gain any rows we inner join so that we only keep the Malaise and Sweep nets filter in alpha_en_data and the samples with black flies/mosquitoes in KBIMP_combined
  filter(!Month %in% c("Jun", "Sep")) %>% 
  group_by(Sector, Sample, Month, Year, Family) %>%
  summarise(Speciessum = n_distinct(Species, na.rm = TRUE), 
            .groups = "drop") %>% #you should go back to double check the issue here to make sure its resolved
  ungroup() %>%
  group_by(Sector, Sample, Year, Month) %>%  
  complete(
    Family = families,
    fill = list(Speciessum = 0)) %>%
  ungroup() 

#number of samples in each group 

num <- speciesrich_CBAYvs_KGLTK %>%
  group_by(Month, Family) %>%
  summarise(unique_samples = n_distinct(Sample), .groups = "drop")

num 

rm(num)

#exploring data for normality and equality of variance

leveneTest(Speciessum ~ Sector, data = speciesrich_CBAYvs_KGLTK) #p-value lesss than 0.05 - variances not equal 
shapiro.test(speciesrich_CBAYvs_KGLTK$Speciessum) #p-value less than 0.05 - data not normal

#running a non-parametric test because the data did not meet the above requirements to determine if species richness varies across sector and month while considering year as an additional factor 

SR_model_bothfam <- art(Speciessum ~ factor(Sector) * factor(Month)  + Error(factor(Year)),
         data = speciesrich_CBAYvs_KGLTK)
anova(SR_model_bothfam)

#separating based on the family 

speciesrich_CBAYvs_KGLTK_bf <- speciesrich_CBAYvs_KGLTK %>%
  filter(Family == "Simuliidae")

SR_model_bf <- art(Speciessum ~ factor(Sector) * factor(Month)  + Error(factor(Year)),
         data = speciesrich_CBAYvs_KGLTK_bf)
anova(SR_model_bf)

#just mosquitoes

speciesrich_CBAYvs_KGLTK_mos <- speciesrich_CBAYvs_KGLTK %>%
  filter(Family == "Culicidae")

SR_model_mos <- art(Speciessum ~ factor(Sector) * factor(Month)  + Error(factor(Year)),
         data = speciesrich_CBAYvs_KGLTK_mos)
anova(SR_model_mos)

#creating data set of total species per month per family for plot

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
  
  facet_wrap(~Family,  
             labeller = labeller(Family = c(
               "Simuliidae" = "Black flies", 
               "Culicidae" = "Mosquitoes")))

ggsave("plots/lineSRtotalSRcbayvskug.png", speciesrichplot , width = 10, height = 4, dpi = 300, bg = "transparent")

#### PART 7 - Multidimentional 2024 and 2025 ----

##### Just mosquitoes -----

KBIMP_speciesmatrix_mosquitoes <- KBIMP_combined %>%
  filter(Family == "Culicidae") %>%
  select(Species, Sample) %>%
  dplyr::count(Sample, Species, name = "presence") %>% 
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

multi_en_data_mos <- meta_data %>%
  filter(Sample %in% rownames(KBIMP_speciesmatrix_mosquitoes)) %>%
  distinct(Sample, .keep_all = TRUE) %>%
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
env_fit <- envfit(NMDS, multi_en_data_mos, permutations = 999, na.rm = TRUE)


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
  distinct(Sample, .keep_all = TRUE) %>%
  filter(SamplingMethod %in% c("Malaise Trap","Sweep Net")) 

KBIMP_speciesmatrix_blackflies <- KBIMP_speciesmatrix_blackflies %>%
  filter(rownames(KBIMP_speciesmatrix_blackflies) %in% multi_en_data_bf$Sample) 

# Order species matrix

KBIMP_speciesmatrix_blackflies <- KBIMP_speciesmatrix_blackflies[order(rownames(KBIMP_speciesmatrix_blackflies)), ]

# Order environmental data

multi_en_data_bf <- multi_en_data_bf %>% arrange(Sample) 

# Bray-Curtis dissimilarity

bray_dist <- vegdist(KBIMP_speciesmatrix_blackflies, method = "bray")

#PERMANOVA

adonis2(bray_dist ~ Sector + Year + SamplingMethod + Month,
        data = multi_en_data_bf, permutations = 999, by="margin")

#multidimensional scaling 

NMDS <- metaMDS(KBIMP_speciesmatrix_blackflies, k=4)

#looking at stress number and plot 

NMDS$stress

plot(
  NMDS$diss,
  NMDS$dist)

#anosim for Sector

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Year

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$Year ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Sampling Method

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$SamplingMethod ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Month

anosim(KBIMP_speciesmatrix_blackflies, multi_en_data_bf$Month ,
       permutations = 999,distance = "bray", strata = NULL)

# Fit environmental vector for CollectionWeek
env_fit <- envfit(NMDS, multi_en_data_bf, permutations = 999, na.rm = TRUE)

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
  
  
  scale_colour_manual(aesthetics = c("colour", "fill"), values = c("CBAY" = "#000099", "KGLTK" =  "#FFC000"), 
                      breaks = c("CBAY", "KGLTK"), labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n", "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"),
                      name = "Region") +
  
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 9, colour = "black"), legend.position = "top")


nmdsplotbf


ggsave("plots/nmdsplotbf2024and202all5.png", nmdsplotbf, width = 6, height = 6, dpi = 300)

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
  distinct(Sample, .keep_all = TRUE) %>%
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

#PERMANOVA

adonis2(bray_dist ~ Sector + Year + SamplingMethod + Month,
        data = multi_en_data_bf, permutations = 999, by="margin")

#Multidimensional scaling 

NMDS <- metaMDS(KBIMP_speciesmatrix_blackflies2, k=3)

#looking at stress number and plot

NMDS$stress

plot(
  NMDS$diss,
  NMDS$dist)

#anosim for Sector

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$Sector ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Year

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$Year ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Sampling Method

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$SamplingMethod ,
       permutations = 999,distance = "bray", strata = NULL)

#anosim for Month

anosim(KBIMP_speciesmatrix_blackflies2, multi_en_data_bf$Month ,
       permutations = 999,distance = "bray", strata = NULL)

# Fit environmental vector

env_fit <- envfit(NMDS, multi_en_data_bf, permutations = 999, 
                  na.rm = TRUE)

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
  distinct(Sample, .keep_all = TRUE) %>%
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

#### PART 8 - Beta diversity analysis  ----

##### comparing the two hamlets #####

#prepare PA matrix

species_matrix_communties <- KBIMP_speciesmatrix_everything %>%
  rownames_to_column(var = "Sample") %>%
  full_join(multi_en_data, join_by(Sample == Sample)) %>%
  select(-Year, -SamplingMethod, -Month, -Lon, -Lat) %>%
  group_by(Sector) %>%
  summarise(across(-Sample, max), .groups = "drop") %>%
  column_to_rownames(var = "Sector")

#run beatdiverity analysis 

betacommunties <- beta.pair(species_matrix_communties, index.family = "sorensen")  

print(betacommunties$beta.sim)
print(betacommunties$beta.sne)
print(betacommunties$beta.sor)

##### temporal beta diversity analysis - Kugluktuk #####

#preparing data 

species_matrix_communties_month <- KBIMP_speciesmatrix_everything %>%
  rownames_to_column(var = "Sample") %>%
  full_join(multi_en_data, join_by(Sample == Sample)) %>%
  filter(Sector == "KGLTK") %>%
  select(-Year, -SamplingMethod, -Sector, -Lat, -Lon) %>%
  group_by(Month) %>%
  summarise(across(-Sample, max), .groups = "drop") %>%
  column_to_rownames(var = "Month")

#run beta diverity analysis 

betacommunties <- beta.pair(species_matrix_communties_month, index.family = "sorensen")  

print(betacommunties$beta.sim)
print(betacommunties$beta.sne)
print(betacommunties$beta.sor)

##### Temporal beta diversity analysis - Cambridge Bay #####

#preparing data 

species_matrix_communties <- KBIMP_speciesmatrix_everything %>%
  rownames_to_column(var = "Sample") %>%
  full_join(multi_en_data, join_by(Sample == Sample)) %>%
  filter(Sector == "CBAY") %>%
  select(-Year, -SamplingMethod, -Sector, -Lon, -Lat) %>%
  group_by(Month) %>%
  summarise(across(-Sample, max), .groups = "drop") %>%
  column_to_rownames(var = "Month")

#run beatdiverity analysis 

betacommunties <- beta.pair(species_matrix_communties, index.family = "sorensen")  

print(betacommunties$beta.sim)
print(betacommunties$beta.sne)
print(betacommunties$beta.sor)

#### PART 9 - Comparing the number of Vectors ----

##### preparing data #####

alpha_en_data <- meta_data %>%
  filter(SamplingMethod %in% c("Malaise Trap", "Sweep Net")) 

families <- c("Culicidae", "Simuliidae") 

vectorrich_CBAYvs_KGLTK <- KBIMP_combined %>%
  inner_join(alpha_en_data) %>% 
  filter(Species %in% vector_change$Species) %>%
  select(-Sequence) %>%
  group_by(Sector, Sample, Year, Month) %>%
  summarise(Speciessum = n_distinct(Species, na.rm = TRUE), 
            .groups = "drop")%>%#you should go back to double check the issue here to make sure its resolved
  ungroup() %>%
  group_by(Sector, Sample, Year, Month) %>%
  complete(
    fill = list(Speciessum = 0)) %>%
  ungroup() 

##### statistical analysis #####

#exploring data for normality and equality of variance

leveneTest(Speciessum ~ Sector, data = vectorrich_CBAYvs_KGLTK ) #p-value lesss than 0.05 - variances not equal 
shapiro.test(vectorrich_CBAYvs_KGLTK$Speciessum) #p-value less than 0.05 - data not normal

#running a non-parametric test because the data did not meet the above requirements to determine if species richness varies across sector and month while considering year as an additional factor 

SR_model_bothfam <- art(Speciessum ~ factor(Sector) + Error(factor(Year)),
                        data = vectorrich_CBAYvs_KGLTK )
anova(SR_model_bothfam)

##### creating a table of the vector species #####

Vectortable <-  vector_change %>%
  select(-Family) %>%
  gt() %>%
  tab_header(title = md("")) %>%
  tab_footnote(footnote = md("")) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(Species, `Vector to`))) %>%
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
            locations = list(cells_column_labels(),
                             cells_row_groups())) %>%
  tab_style(style = cell_borders(sides = c("bottom"),
                                 color = "black", weight = px(2)),
            locations = cells_title()) %>%
  tab_style(style = cell_borders(sides = c("top"),
                                 color = "black", weight = px(2)),
            locations = cells_footnotes()) %>%
  tab_style(style = cell_borders(
    sides = c("bottom"),
    color = "black", weight = px(2)),
    locations = cells_column_labels()) %>%
  tab_options(data_row.padding = px(5)) 

gtsave(data = Vectortable, 
       filename = "plots/Vectortable.png")

##### iNEXT for confiedence intervals #####

#getting the confidence intervals for vectors 

iNEXT_vectors <- KBIMP_combined %>%
  select(Sample, Species) %>%
  mutate(Species = if_else(Species == "Aedes nigripes/impiger", 
                           "Aedes nigripes/Aedes impiger", Species)) %>%
  mutate(Species = str_replace(Species, "Simulium arcticum complex sp", 
                               "Simulium arcticum complex")) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  filter(Species %in% vector_change$Species) %>%
  add_row(Species = "sampling_extent",
          CBAY = 209,
          KGLTK = 114, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext_vectors <- iNEXT(iNEXT_vectors, q=0, datatype="incidence_freq")

em.inext_vectors$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_vectors <- as.data.frame(em.inext_vectors$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(type = "Vectors")

em.inext_vectors$DataInfo

#getting the confidence intervals for non-vectors 

iNEXT_nonvectors <- KBIMP_combined %>%
  select(Sample, Species) %>%
  mutate(Species = if_else(Species == "Aedes nigripes/impiger", "Aedes nigripes/Aedes impiger", Species)) %>%
  mutate(Species = str_replace(Species, "Simulium arcticum complex sp", "Simulium arcticum complex")) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  filter(!Species %in% vector_change$Species) %>%
  add_row(Species = "sampling_extent",
          CBAY = 209,
          KGLTK = 114, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext_nonvectors <- iNEXT(iNEXT_nonvectors, q=0, datatype="incidence_freq")

em.inext_nonvectors$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_nonvectors <- as.data.frame(em.inext_nonvectors$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(type = "Non-vectors")

em.inext_nonvectors$DataInfo

#combining the data set for vectors, non-vectors, and total

vector2024 <- bind_rows(confidenceinterval_nonvectors, 
                        confidenceinterval_total, 
                        confidenceinterval_vectors)

vector2024 <- vector2024 %>%
  mutate(Year = "2024-2025")

##### Figure production for vectors #####

vectorchange <- ggplot() +
  
  geom_point(data = vector2024,
             aes(x = Sector, y = qD, colour = type, group = type), 
             size = 3) +
  
  geom_errorbar(data = vector2024, aes(x = Sector, ymin = qD.LCL, 
                                       ymax = qD.UCL, colour = type),
                width = 0.15) +
  
  scale_colour_manual(
    values = c("Vectors" = "#9D5C63", "Non-vectors" =  "#78BC61", "Total" = "black"))  +
  
  scale_x_discrete(
    breaks = c("CBAY", "KGLTK"),
    labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n",
                "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"), name = "Region") +
  
  theme_bw(base_size = 15) +
  xlab("Sampling Years") +
  ylab("Total Species Richness") +
  
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14)) 

ggsave("plots/vector.png", vectorchange , width = 8, height = 4, dpi = 300, bg = "transparent")


#### PART 10 - Comparing the number of blood feeding species  ----

##### preparing data #####

families <- c("Culicidae", "Simuliidae") 

speciesrich_CBAYvs_KGLTK_nb <- KBIMP_combined %>%
  inner_join(alpha_en_data) %>% #I think I fixed the problem here so we do not gain any rows we inner join so that we only keep the Malaise and Sweep nets filter in alpha_en_data and the samples with black flies/mosquitoes in KBIMP_combined
  filter(!Month %in% c("Jun", "Sep")) %>% 
  filter(!Species %in% nonbiting_species$Species) %>%
  group_by(Sector, Sample, Month, Year, Family) %>%
  summarise(Speciessum = n_distinct(Species, na.rm = TRUE), 
            .groups = "drop") %>% #you should go back to double check the issue here to make sure its resolved
  ungroup() %>%
  group_by(Sector, Sample, Year, Month) %>%  
  complete(
    Family = families,
    fill = list(Speciessum = 0)) %>%
  ungroup() 

#number of samples in each group 

num_nb <- speciesrich_CBAYvs_KGLTK %>%
  group_by(Month, Family) %>%
  summarise(unique_samples = n_distinct(Sample), .groups = "drop")

num_nb 

rm(num_nb)

##### statistical analysis #####

#exploring data for normality and equality of variance

leveneTest(Speciessum ~ Sector, data = speciesrich_CBAYvs_KGLTK_nb) #p-value lesss than 0.05 - variances not equal 
shapiro.test(speciesrich_CBAYvs_KGLTK_nb$Speciessum) #p-value less than 0.05 - data not normal

#running a non-parametric test because the data did not meet the above requirements to determine if species richness varies across sector and month while considering year as an additional factor 

#separating based on the family 

speciesrich_CBAYvs_KGLTK_bf_nb <- speciesrich_CBAYvs_KGLTK_nb %>%
  filter(Family == "Simuliidae")

SR_model_bf_nb <- art(Speciessum ~ factor(Sector)  + Error(factor(Year)),
                   data = speciesrich_CBAYvs_KGLTK_bf_nb)
anova(SR_model_bf_nb)

##### iNEXT for confidence intervals #####

#getting the confidence intervals for non-biting

iNEXT_nonbiting <- KBIMP_combined %>%
  select(Sample, Species) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  filter(Species %in% nonbiting_species$Species) %>%
  add_row(Species = "sampling_extent",
          CBAY = 209,
          KGLTK = 114, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext_nonbiting <- iNEXT(iNEXT_nonbiting, q=0, datatype="incidence_freq")

em.inext_nonbiting$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_nonbiting <- as.data.frame(em.inext_nonbiting$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(type = "Non-biting")

em.inext_nonbiting$DataInfo

#getting the confidence intervals for non-vectors 

iNEXT_biting <- KBIMP_combined %>%
  filter(Family == "Simuliidae") %>%
  select(Sample, Species) %>%
  mutate(Species = str_replace(Species, "Aedes punctor/Aedes hexodontus", "Aedes hexodontus")) %>%
  mutate(Species = if_else(Species == "Aedes nigripes/impiger", "Aedes nigripes/Aedes impiger", Species)) %>%
  mutate(Species = str_replace(Species, "Simulium arcticum complex sp", "Simulium arcticum complex")) %>%
  mutate(region = str_extract(Sample, "^[A-Za-z]+")) %>%
  dplyr::count(region, Sample, Species, name = "Abundance") %>% 
  mutate(Abundance = ifelse(Abundance > 0, 1, 0)) %>%
  dplyr::count(region, Species, name = "Incidence") %>%
  pivot_wider(names_from = region, values_from = Incidence) %>%
  mutate(CBAY = replace_na(CBAY, 0), KGLTK = replace_na(KGLTK, 0)) %>%
  filter(!Species %in% nonbiting_species$Species) %>%
  add_row(Species = "sampling_extent",
          CBAY = 209,
          KGLTK = 114, .before = 1) %>%
  column_to_rownames(var = "Species")

em.inext_biting <- iNEXT(iNEXT_biting, q=0, datatype="incidence_freq")

em.inext_biting$iNextEst

#isolating the data on confidence intervals for the observed data set 

confidenceinterval_biting <- as.data.frame(em.inext_biting$iNextEst$size_based) %>%
  filter(Method == "Observed") %>%
  select(Assemblage, qD, qD.LCL, qD.UCL) %>%
  mutate(x = "total") %>%
  rename_with(~"Sector", contains("Assemblage")) %>%
  mutate(type = "Biting")

em.inext_biting$DataInfo

#combining the biting, non-biting, and total data sets 

biting2024 <- bind_rows(confidenceinterval_biting, 
                        confidenceinterval_total, 
                        confidenceinterval_nonbiting)

biting2024 <- biting2024 %>%
  mutate(Year = "2024-2025") %>%
  mutate(type = replace_na(type, "Total"))

##### figure for biting vs non-biting species #####

bitingchange <- ggplot() +
  
  geom_point(data = biting2024,
             aes(x = Sector, y = qD, colour = type, group = type), size = 3) +
  
  geom_errorbar(data = biting2024, aes(x = Sector, ymin = qD.LCL, ymax = qD.UCL, colour = type),
                width = 0.15) +
  
  scale_colour_manual(
    values = c("Biting" = "#9D5C63", "Non-biting" =  "#78BC61", "Total" = "black"))  +
  
  scale_x_discrete(
    breaks = c("CBAY", "KGLTK"),
    labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n",
                "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"), name = "Region") +
  
  theme_bw(base_size = 15) +
  xlab("Sampling Years") +
  ylab("Total Species Richness") +
  
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14)) 

ggsave("plots/bitingchange.png", bitingchange , width = 8, height = 4, dpi = 300, bg = "transparent")


##### combining vector and biting plot #####

combinedchangeplot <- (bitingchange + vectorchange) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = list("A"))

png("plots/changeplot.png", width =3000, height = 1500, res = 300)

print(combinedchangeplot)

dev.off()

#### PART 11 - Map of regions and samplng locations  ----

##### Map for Kugluktuk sample locations #####

#getting lat/lon from our metadata set 

kugluktuk_mapdata <- meta_data %>%
  filter(Sector == "KGLTK") %>%
  mutate(Lat = clean_coords(Lat)) %>%
  mutate(Lon = clean_coords(Lon)) %>%
  mutate(Lat = gsub("^[^ ]* ", "", Lat), 
         Lat = as.numeric(Lat), Lat = if_else(Sector == "KGLTK", 67 + (Lat / 60), Lat), 
         Lon = gsub("^[^ ]* ", "", Lon), Lon = as.numeric(Lon),
         Lon = if_else(Sector == "KGLTK", -(115 + (Lon / 60)), Lon)) %>%
  filter(!is.na(Lat), !is.na(Lon)) %>%
  mutate(Year = as.character(Year)) %>%
  filter(!Sample == "KGLTK0271")

#reformatting this for maps

Coord_images<- kugluktuk_mapdata %>%
  group_by(Lon, Lat)%>%
  summarise(coord= n())%>%
  filter(!is.na(Lat))%>%
  filter(!is.na(Lon))

#getting the world map dataset from maps and creating the base map

world_map <- map_data("world")

map_data_cr <- map_data('world')[map_data('world')$region == "Nunavut",]
p <- ggplot() + coord_fixed(1.3, xlim = c(-107, -105), ylim = c(65, 69))


base_world_messy <- p + 
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
               colour="Gainsboro", fill="Gainsboro")+
  geom_polygon(data = map_data_cr, aes(x=long, y=lat, group = group),
               colour = 'red', fill = 'pink')

cleanup <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = 'white', colour = 'white'), 
                 axis.line = element_line(colour = "white"), legend.position="none",
                 axis.ticks=element_blank(), axis.text.x=element_blank(),
                 axis.text.y=element_blank())

base_world <- base_world_messy + cleanup
base_world+
  geom_point(data=Coord_images, aes(x= Lon, y= Lat), colour="Dim Gray", fill="Dim Gray", size=3)

#getting our lat/lon into map format 

ll_means <- sapply(kugluktuk_mapdata[5:6], mean)
sq_map2 <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 11)

#plotting sample points on map 

map_kug <- ggmap(sq_map2) + 
  geom_point(data = kugluktuk_mapdata , aes(x = Lon, y = Lat, colour = Sector), size = 2, shape = 19) +
  scale_color_manual(values = c("#FFC000" )) +  
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis label
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"))   # Remove y-axis label

map_kug

ggsave("plots/map_kugluktuk.png", map_kug, width = 4, height = 3, dpi = 300)

##### Map for Cambridge Bay sample locations #####

#getting lat/lon from our metadata set 

cbay_mapdata <- meta_data  %>%
  filter(Sector == "CBAY") %>%
  filter(!is.na(Lat), !is.na(Lon)) %>%
  mutate(Lon = as.numeric(Lon)) %>%
  mutate(Lat = as.numeric(Lat))

#reformatting this for maps

Coord_images<- cbay_mapdata %>%
  group_by(Lon, Lat)%>%
  summarise(coord= n())

#creating world base map 

map_data_cr <- map_data('world')[map_data('world')$region == "Nunavut",]
p <- ggplot() + coord_fixed(1.3, xlim = c(-106, -104), ylim = c(67, 70))

base_world_messy <- p + 
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
               colour="Gainsboro", fill="Gainsboro")+
  geom_polygon(data = map_data_cr, aes(x=long, y=lat, group = group),
               colour = 'red', fill = 'pink')

cleanup <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = 'white', colour = 'white'), 
                 axis.line = element_line(colour = "white"), legend.position="none",
                 axis.ticks=element_blank(), axis.text.x=element_blank(),
                 axis.text.y=element_blank())

base_world <- base_world_messy + cleanup + 
  geom_point(data=Coord_images, aes(x= Lon, y= Lat), colour="Dim Gray", fill="Dim Gray", size=3)

#transforming our data into the correct format 

ll_means <- sapply(cbay_mapdata[5:6], mean)
sq_map2 <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 10)

#creating our map of sample points

mapcbay <- ggmap(sq_map2) + 
  geom_point(data = cbay_mapdata , aes(x = Lon, y = Lat, colour = Sector), size = 2, shape = 19) +
  scale_color_manual(values = c("#000099")) +  
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis label
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"))   # Remove y-axis label

mapcbay

ggsave("plots/map_cbay.png", mapcbay, width = 4, height = 3, dpi = 300)


##### Nunavut locations map #####

#getting lat/lon coordinates for cbay and kgltk

both_mapdata <- meta_data  %>%
  filter(Sample %in% c("KGLTK0056", "CBAY0040")) %>%
  mutate(Lat = clean_coords(Lat)) %>%
  mutate(Lon = clean_coords(Lon)) %>%
  mutate(Lat = gsub("^[^ ]* ", "", Lat), 
         Lat = as.numeric(Lat), Lat = if_else(Sector == "KGLTK", 67 + (Lat / 60), Lat), 
         Lon = gsub("^[^ ]* ", "", Lon), Lon = as.numeric(Lon),
         Lon = if_else(Sector == "KGLTK", -(115 + (Lon / 60)), Lon)) %>%
  mutate(Lon = as.numeric(Lon)) %>%
  mutate(Lat = as.numeric(Lat))

#direct shapefile URL for 2021 Census Divisions for the Kitikmeot region and Nunavut

st_can_url <- "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lcd_000b21a_e.zip"

#download and unzip automatically

temp_zip <- tempfile(fileext = ".zip")
temp_dir <- tempdir()
download.file(st_can_url, temp_zip, quiet = TRUE)
unzip(temp_zip, exdir = temp_dir)

#load boundary and filter for Kitikmeot (CDUID 6208)

kitikmeot_sf <- st_read(file.path(temp_dir, "lcd_000b21a_e.shp"), quiet = TRUE) %>% 
  filter(CDUID == "6208") %>% 
  st_transform(4326)

#load boundary and filter for Nunavut (PRUID 62)

nunavut_sf <- st_read(file.path(temp_dir, "lcd_000b21a_e.shp"), quiet = TRUE) %>% 
  filter(PRUID == "62") %>% 
  st_transform(4326)

# Your points as sf
points_sf <- st_as_sf(both_mapdata, coords = c("Lon", "Lat"), crs = 4326)


#create the Nunavut map

nunavut_map <- ggplot() +
  # Draw Nunavut boundary
  geom_sf(data = nunavut_sf, fill = "#e9ecef", color = "#495057", linewidth = 0.3) +
  
  geom_sf(data = kitikmeot_sf, fill = "#A2CD5A", color = "#1b4965", linewidth = 0.4, alpha = 0.7) +
  
  geom_sf(data = points_sf, aes(color = Sector), size = 2, alpha = 0.9) +

  
  # Set an Arctic-appropriate projection (EPSG:3573)
  coord_sf(crs = 3573) +
  
  scale_colour_manual(
    values = c("CBAY" = "#000099", "KGLTK" = "#FFC000"),
    breaks = c("CBAY", "KGLTK"),
    labels = c(
      "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n",
      "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"
    ), name = "Region") +
  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    color = NA), 
    panel.grid = element_line(color = "white", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_blank(),
    axis.text = element_blank())

# Save and export map
png("plots/nunavut_sampling_map2.png", width = 1000, height = 800, res = 300)
print(nunavut_map)
dev.off()

##### Circumpolar map #####

Coord_images<- both_mapdata %>%
  group_by(Lon, Lat)%>%
  summarise(coord= n())

# World map as sf
world_sf <- st_as_sf(map_data("world"), coords = c("long", "lat"), crs = 4326, agr = "constant") %>%
  group_by(group) %>%
  summarise(do_union = FALSE) %>%
  st_cast("POLYGON")

#creating arctic circle "point" to put on map

arctic_circle <- data.frame(
  Lon = seq(-180, 180, length.out = 1000),
  Lat = 66.56)

arctic_circle_sf <- st_as_sf(arctic_circle, coords = c("Lon", "Lat"), crs = 4326) %>%
  summarise() %>%
  st_cast("LINESTRING")

#formatting map to be of the arctic circle 

polar_crs <- 3995  # Arctic polar stereographic

#get clean world polygons (better than map_data)

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

#keep only Arctic region (clip at 66.56°N)

arctic_boundary <- st_sfc(
  st_polygon(list(cbind(
    c(-180, 180, 180, -180, -180),
    c(66.56, 66.56, 90, 90, 66.56)))),
  crs = 4326)

world_arctic <- st_intersection(world_sf, arctic_boundary)

world_arctic <- st_crop(
  world_sf,
  xmin = -180, xmax = 180,
  ymin = 60, ymax = 90)

#clean world polygons

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

#convert to lat/lon explicitly

world_sf <- st_transform(world_sf, 4326)

#keep only features that extend into Arctic

world_arctic <- world_sf %>%
  filter(st_bbox(.)[["ymax"]] >= 66.56)

#points

points_sf <- st_as_sf(both_mapdata, coords = c("Lon", "Lat"), crs = 4326)

#arctic Circle line

arctic_circle <- st_sfc(
  st_linestring(cbind(seq(-180, 180, length.out = 1000), rep(66.56, 1000))),
  crs = 4326)

#plotting the Circumpolar Map

circumpolarmap <- ggplot() +
  #base Arctic countries
  geom_sf(data = world_arctic, fill = "grey90", color = "grey40", linewidth = 0.2) +
  
  #highlight Kitikmeot Region
  geom_sf(data = kitikmeot_sf, fill = "#A2CD5A", color = "#1b4965", linewidth = 0.4, alpha = 0.7) +
  
  #arctic Circle line (Fixed)
  geom_sf(data = arctic_circle, color = "black", linewidth = 0.8, linetype = "solid") +
  
  #CBAY and KGLTK
  geom_sf(data = points_sf, aes(color = Sector), size = 3, alpha = 0.9) +
  
  #polar projection limits
  coord_sf(crs = 3995, ylim = c(-3e6, 3e6), xlim = c(-3e6, 3e6)) +
  
  scale_colour_manual(
    values = c("CBAY" = "#000099", "KGLTK" = "#FFC000"),
    breaks = c("CBAY", "KGLTK"),
    labels = c(
      "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n",
      "KGLTK" = "Kugluktuk\n(Qurluqtuq)\n"), 
    name = "Region") +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank())

png("plots/circumpolarmap.png", width = 1400, height = 1400, res = 300)

print(circumpolarmap)

dev.off()

#### PART 12 - Temperature data figure ----

##### Cambridge Bay #####

#loading temperature data from the list of files

temp_cbay <- list.files(
  path = "raw-data2/temp_data_cbay",
  pattern = "\\.csv$", full.names = TRUE)

#combining all these files into one

temp_cbay_data <- do.call(rbind, lapply(temp_cbay, read.csv))

#finding the summer temp (July and August)  

temp_cbay_avgmonthly <- temp_cbay_data %>%
  select(Date.Time, Year, Month, Max.Temp...C., Min.Temp...C., Mean.Temp...C.) %>%
  filter(Month %in% c( "7", "8")) %>%
  group_by(Year) %>%
  summarise(avg_monthly_temp = mean(Mean.Temp...C., na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Sector = "CBAY")


ggplot(temp_cbay_avgmonthly, aes(x = Year, y = avg_monthly_temp)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "red")

#creating the model line of temp increase 

model <- lm(avg_monthly_temp ~ Year, data = temp_cbay_avgmonthly)
temp_cbay_avgmonthly$trend <- predict(model)

###### Kugluktuk #####

#loading temperature data from the list of files

temp_kug <- list.files(
  path = "raw-data2/temp_data_kug",
  pattern = "\\.csv$", full.names = TRUE)

#combining all these files into one

temp_kug_data <- do.call(rbind, lapply(temp_kug, read.csv))

#finding the summer temp (July and August)  

temp_kug_avgmonthly <- temp_kug_data %>%
  select(Date.Time, Year, Month, Max.Temp...C., Min.Temp...C., Mean.Temp...C.) %>%
  filter(Month %in% c( "7", "8")) %>%
  group_by(Year) %>%
  summarise(avg_monthly_temp = mean(Mean.Temp...C., na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Sector = "KUG")


ggplot(temp_kug_avgmonthly, aes(x = Year, y = avg_monthly_temp)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "red")

#creating the model line of temp increase '

model <- lm(avg_monthly_temp ~ Year, data = temp_kug_avgmonthly)
temp_kug_avgmonthly$trend <- predict(model)

##### combined plot #####

#combine the data sets for both locations of average temp

temp_avgmonthly <- bind_rows(temp_cbay_avgmonthly, temp_kug_avgmonthly)

#plotting the data 

increasingtemp <- ggplot(temp_avgmonthly, aes(x = Year, y = avg_monthly_temp, 
                                              colour = Sector, group = Sector)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(
    values = c("CBAY" = "#000099", "KUG" =  "#FFC000"),
    breaks = c("CBAY", "KUG"),
    labels = c( "CBAY"  = "Cambridge Bay\n(Iqaluktuuttiaq)\n",
                "KUG" = "Kugluktuk\n(Qurluqtuq)\n"), name = "Region") +
  
  theme_bw(base_size = 15) +
  ylab("Average Summer Temperature") +
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14)) 


ggsave("plots/increasingtemps.png", increasingtemp , width = 8, height = 4, dpi = 300, bg = "transparent")


#### PART 13 - Creating range maps for vector species ----

##### preparing ecozone data #####

world <- ne_download(scale = "large", type = "land", category = "physical", returnclass = "sf")

canada <- ne_countries(country = "Canada", returnclass = "sf")

url <- "https://services.arcgis.com/lGOekm0RsNxYnT3j/ArcGIS/rest/services/National_ecological_framework_of_Canada_ecozones/FeatureServer/0/query?where=1=1&outFields=*&outSR=4326&f=geojson"

ecozones <- st_read(url)


land <- ne_download(
  scale = "medium",
  type = "land",
  category = "physical",
  returnclass = "sf")

canada <- ne_countries(country = "Canada", returnclass = "sf")

# Clip land to Canada
canada_land <- st_intersection(land, canada)


canada_land <- st_transform(canada_land, st_crs(ecozones))

ecozones_land <- st_intersection(ecozones, canada_land)

ecozones2 <- ecozones_land %>%
  mutate(ECOZONE_NAME_EN = case_when(ECOZONE_NAME_EN 
                                     %in% c("Taiga Shield",
                                            "Taiga Plains",
                                            "Taiga Cordillera")
                                     ~ "Taiga" , 
                                     ECOZONE_NAME_EN 
                                     %in% c("Southern Arctic",
                                            "Arctic Cordillera",
                                            "Northern Arctic") 
                                     ~ "Arctic",
                                     ECOZONE_NAME_EN 
                                     %in% c("Boreal Cordillera",
                                            "Boreal Plains",
                                            "Boreal Shield") 
                                     ~ "Boreal",
                                     ECOZONE_NAME_EN 
                                     %in% c("Atlantic Maritime",
                                            "Pacific Maritime") 
                                     ~ "Maritime", 
                                     TRUE ~ ECOZONE_NAME_EN))

taiga <- ecozones_land %>%
  filter(ECOZONE_NAME_EN %in% c(
    "Taiga Shield",
    "Taiga Plains",
    "Taiga Cordillera")) %>%
  summarise() 

arctic <- ecozones_land %>%
  filter(ECOZONE_NAME_EN %in% c(
    "Southern Arctic", 
    "Arctic Cordillera"
  )) %>%
  summarise()

arctic <- st_boundary(arctic)
taiga <- st_boundary(taiga)
taiga <- st_make_valid(taiga)
arctic <- st_make_valid(arctic)

treeline <- st_intersection(
  st_union(st_buffer(taiga, 3000)),
  st_union(st_buffer(arctic, 3000))
) %>%
  st_boundary()


my_colors <- paletteer_d("ggthemes::Tableau_20")

##### organizing my data ##### 

recordsmap_mydata <- KBIMP_combined %>%
  left_join(meta_data) %>%
  mutate(decimalLatitude = round(Lat, 0)) %>%
  mutate(decimalLongitude = round(Lon, 0)) %>%
  select(Species, decimalLongitude, decimalLatitude, Year) %>%
  dplyr::rename(species = Species) %>%
  dplyr::rename(year = Year) %>%
  mutate(data_type = "Our data")

##### Aedes excrucians #####

#find the user key for Insecta (with gbif it is likely easier to download all of the data for insecta in Canada and then filter for the required species groups rather than creating many usage keys)

aedesexcrucians_binf <- records_gbif_mos %>%
  filter(species == "Aedes excrucians") %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)

aedesexcrucians_ourdata <- recordsmap_mydata %>%
  filter(species == "Aedes excrucians") %>%
  distinct() 
  
aedesexcrucians_histdata <- bind_rows(aedesexcrucians_binf, aedesexcrucians_ourdata)

aedesexcrucians_histdata <- aedesexcrucians_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))


aedesexcruciansmap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = aedesexcrucians_histdata , aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Schaefer" = "Schaefer 2012",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source", color = "Time Frame", title = "Aedes excrucians") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))


##### Simulium decorum #####

decorum_binf  <- records_gbif_bf  %>%
  filter(species == "Simulium decorum") %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)

simuliumdecorum_ourdata <- recordsmap_mydata %>%
  filter(species == "Simulium decorum") %>%
  distinct() 

decorum_histdata <- bind_rows(decorum_binf, simuliumdecorum_ourdata)

decorum_histdata <- decorum_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))

decorummap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = decorum_histdata , aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20, "Schaefer" = 18), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Schaefer" = "Schaefer 2012",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source",  color = "Time Frame",title = "Simulium decorum") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))


##### Aedes communis #####

aedescommunis_binf  <- records_gbif_mos  %>%
  filter(species == "Aedes communis") %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)
  
aedescommunis_ourdata <- recordsmap_mydata %>%
    filter(species == "Aedes communis") %>%
    distinct() 

aedescommunis_histdata <- bind_rows(aedescommunis_binf, aedescommunis_ourdata)

aedescommunis_histdata <- aedescommunis_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))

aedescommunismap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = aedescommunis_histdata , aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20, "Schaefer" = 18), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Schaefer" = "Schaefer 2012",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source", color = "Time Frame", title = "Aedes communis") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))


##### Simulium vittatum #####

vittatum_binf  <- records_gbif_bf  %>%
  filter(species == "Simulium vittatum") %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)

vittinum_ourdata <- recordsmap_mydata %>%
    filter(species == "Simulium vittatum") %>%
    distinct() 

vittinum_histdata <- bind_rows(vittatum_binf, vittinum_ourdata)

vittinum_histdata <- vittinum_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))

vittinummap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = vittinum_histdata , aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20, "Schaefer" = 18), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Schaefer" = "Schaefer 2012",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source", color = "Time Frame", title = "Simulium vittatum") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))

##### Simulium norelli #####

noelleri_binf  <- records_gbif_bf  %>%
  filter(species == "Simulium noelleri") %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)

noelleri_ourdata <- recordsmap_mydata %>%
  filter(species == "Simulium noelleri") %>%
  distinct()

noelleri_histdata <- bind_rows(noelleri_binf, noelleri_ourdata)

noelleri_histdata <- noelleri_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))

noellerimap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = noelleri_histdata , aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20, "Schaefer" = 18), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Schaefer" = "Schaefer 2012",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source", color = "Time Frame", title = "Simulium noelleri") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))


##### Aedes hexodontus/ Aedes punctor #####

hexodontus_binf  <- records_gbif_mos  %>%
  filter(species %in% c("Aedes hexodontus","Aedes punctor")) %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)

hexodontus_ourdata <- recordsmap_mydata %>%
  filter(species == "Aedes punctor/Aedes hexodontus") %>%
  distinct()

hexodontus_histdata <- bind_rows(hexodontus_binf, hexodontus_ourdata)

hexodontus_histdata <- hexodontus_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))

hexodontusmap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = hexodontus_histdata , aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20, "Schaefer" = 18), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source", color = "Time Frame", title = "Aedes hexodontus/\nAedes punctor\n") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))


##### Aedes nigripis/ Aedes impiger #####

nigripis_binf  <- records_gbif_mos  %>%
  filter(species %in% c("Aedes nigripes","Aedes impiger")) %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)

nigripis_ourdata <- recordsmap_mydata %>%
  filter(species == "Aedes nigripes/impiger") %>%
  distinct()

nigripis_histdata <- bind_rows(nigripis_binf, nigripis_ourdata)

nigripis_histdata <- nigripis_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))

nigripismap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = nigripis_histdata, aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source", color = "Time Frame", title = "Aedes nigripes/\nAedes impiger\n") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))


##### Simulium norelli #####

venustum_binf  <- records_gbif_bf  %>%
  filter(species == "Simulium venustum") %>%
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) %>%
  mutate(data_type = "GBIF") %>%
  select(decimalLatitude, decimalLongitude, data_type, year)

venustum_ourdata <- recordsmap_mydata %>%
  filter(species == "Simulium venustum complex") %>%
  distinct()

venustum_histdata <- bind_rows(venustum_binf, venustum_ourdata)

venustum_histdata <- venustum_histdata %>%
  mutate(period = ifelse(year <= 2015, "Before 2015", "2015-2026"))

venustummap <- ggplot() +
  geom_sf(data = ecozones2, aes(fill = ECOZONE_NAME_EN)) +
  geom_sf(data = treeline, color = "black", size = 1.2) +
  geom_point(data = venustum_histdata , aes(x = decimalLongitude, y = decimalLatitude, colour = period, shape = data_type), 
             stroke = 1, size = 3) +
  scale_fill_manual(values= my_colors) +
  scale_colour_manual(values = c("2015-2026" = "red", "Before 2015" = "black")) +
  scale_shape_manual(values = c("Our data" = 9, "GBIF" = 20, "Schaefer" = 18), 
                     labels = c( "GBIF"  = "GBIF",
                                 "Schaefer" = "Schaefer 2012",
                                 "Our data" = "Our Data")) +
  
  labs( fill = "Ecozone", shape = "Data source", color = "Time Frame", title = "Simulium venustum") +
  
  theme_bw() +
  theme(title = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 14))



##### combing plots #####

combined_plot <- (aedesexcruciansmap + aedescommunismap + hexodontusmap+ nigripismap + decorummap + vittinummap + noellerimap + venustummap) +
  plot_layout(guides = "collect", ncol = 4) 

png("plots/histmaps.png", width = 4000, height = 2000, res = 300)

print(combined_plot)

dev.off()

#### PART 14 - Positive control figure ----

poscontrolplatemap <- read_csv(file = "raw-data2/poscontrolplatemap.csv")

POSCON_seqs <- read_tsv(file = "raw-data2/Shauna_COINEM_POS_TaxonomicAssignments_DominantContigs.tsv")


POSCON_seqs <- POSCON_seqs %>%
  select(Sample, Species, ReadCount) %>%
  right_join(poscontrolplatemap, join_by(Sample == Well)) %>%
  distinct() %>%
  filter(!is.na(Species)) %>%
  pivot_longer(cols = -c(Sample, Species, ReadCount, "Black fly concentration", "Mosquito concentration"), 
               names_to = "Nematode Species", 
               values_to = "Nem_Concentration") %>%
  pivot_longer(cols = -c(Sample, Species, ReadCount, "Nematode Species", "Nem_Concentration"), 
               names_to = "Insect Species", 
               values_to = "Insect Concentration") %>%
  mutate(across(c(ReadCount, `Insect Concentration`, `Nem_Concentration`), as.numeric)) %>%
  filter(`Nem_Concentration` != 0) %>%
  mutate(`Insect Concentration` = log(`Insect Concentration` + 0.000001),   `Nem_Concentration` = log(`Nem_Concentration` + 0.0000001))



bubbleplot <-ggplot(POSCON_seqs, aes(x = `Insect Concentration`, 
                                     y= `Nem_Concentration` , 
                                     size = ReadCount, 
                                     colour = `Nematode Species`)) +
  geom_point() +
  xlim(-15, 3) +
  ylim(-6, 3) +
  scale_colour_manual(values = c("Onchacerca concentration" = "#9DD1F1", "Seteria concentration" = "#E65F5C")) +
  theme_bw() +
  facet_grid(~ `Nematode Species`) +
  theme(axis.text.x = element_text(angle = 0, size = 10, color = "black"),
        axis.text.y = element_text(angle = 0, size = 10, color = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.background = element_rect(fill = NA, color = NA), legend.background = element_rect(fill = NA, color = NA), 
        strip.background = element_rect(fill = NA, color = NA),
        strip.text  = element_text(face = "bold", size = 12)) 


ggsave("plots/bubbleplot.png", bubbleplot , width = 8, height = 3, dpi = 300, bg = "transparent")








