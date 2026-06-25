# 📊 Project overview 


In this study we used sequences generated from DNA barcodng and metabarcoding on the Oxford Nanopore MinION 
to investigate nematode incidence and biting insect diversity in the Arctic. The data has already been run through 
The Barcode Inference Pipeline (BIP) which was developed at the centre for biodiverity genomics. Here we filter out contmanination 
detected on the negatives controls and identify species which could not be identified with BIP (which uses probabolistic methods) 
using sequence simularity on BOLD. To investigate the successfulness of this method we build phylogentic trees using maximum liklihood 
and 1000 bootstrapping iterations. We then investigated species diversity using a variatey of packages such as iNEXT, vegan, and betapart. 

6 script files, only 1, 2, 5, and 6 are needed to run the analysis (3 and 4 are another method of species assignemnt). 

## 📁 Script Files

### 1. `1_Data_Cleaning.R`
- Contains all data cleaning steps recorded in the lab notebook  
- Includes:
  - Fixes for problematic samples  
  - Handling of samples rerun on new plates  
  - Contamination correction using negative controls  
- Processing is separated by year:
  - **2024 (barcoding)**  
  - **2025 (metabarcoding)** (different control plate setup)

---

### 2. `2_Functions.R`
- Includes all custom functions used across analysis scripts  

---

### 3. `3_Phylogenetic_Analysis_Mediod_2024.R`
- Phylogenetic analysis for **2024 samples**  
- Includes:
  - Tree construction  
  - Species assignment for unresolved sequences  
- Uses the **medoid method**

---

### 4. `4_Phylogenetic_Analysis_Mediod_2025.R`
- Same workflow as above, but for **2025 samples**  
- Uses the **medoid method**

---

### 5. `5_Phylogenetic_Analysis_UniqueSequences.R`
- Combined phylogenetic analysis for **2024 + 2025**
- Uses **unique sequences**
- Includes:
  - BOLD database sequence import  
  - Comparison of phylogenetic assignment methods (Scripts 3 & 4)

---

### 6. `6_Ecological_Analysis.R`
- Contains ecological analyses and most figure generation:
  - iNEXT analyses  
  - Comparison with **2012 dataset**  
  - Species richness calculations  
  - NMDS  
  - GBIF map generation  
  - Sample map production  
  - Temperature comparison plots  
- Clearly divided into labeled sections  

---

## 📂 Data

### 🔹 Raw Data (`Raw-data2`)
Contains all primary data used in analysis.

#### 2024 Sequence Data
- `KGLKTK_2024_OTUDetails.tsv`  
- `CBAY2024_AllPlates_OTUDetails.tsv`  
  - Includes `problemsamples.csv` (lab-noted issues)  
- `Shauna_CBAY2024_Plate3_OTUDetails.tsv`

#### 2025 Sequence Data
- `KBIMP2025_insectCOI_OTUDetails.tsv`  
- Control files:
  - `extractiondata2025.csv`  
  - `PCRcontrol_data.csv`

#### Outgroup Data
- `Outgroup.csv`

#### Metadata Files
- `CBAY2025_metadata.csv`  
- `KGLTK2025_metadata.csv`  
- `condencedsites.csv`  
- `KBIMP2024_specimendata.csv`  
- `KBIMP_meta_sitenamesfixed.csv`  
- `KBIMP2024_abundence.csv`  
- `vector_change.csv`

#### Historical Dataset
- `schafer_2012.csv`

---

### 🔹 Processed Data
Generated during analysis and reused in later steps.

#### Filtered Sequence Files
- `KBIMP2024_filteredCOI.tsv`  
- `KBIMP2025_filteredCOI.tsv`

#### Unique Sequence Files (for BOLD searches)
- `uniquemosseqforbold.fasta`  
- `unique_simseq_forbold.fasta`  
- `unique_bfseq_forbold.fasta`

#### BOLD Database Outputs
- `BOLDaedescombined.csv`  
- `BOLDID_simuulidae.csv`  
- `BOLDID_notsimulidae.csv`  
- `BOLDIDspecies.csv`  
  - Contains **BIN codes and corresponding species**

#### Final Species Dataset
- `KBIMP_updatedspecies.tsv`

---

## 📂 Required packages 

The analysis scripts require R and the following packages (versions used in the manuscript are noted):

stringr
tidyverse
readr
viridis
ggplot2 
Biostrings 
ape 
muscle
phangorn    
ggtree    
seqinr
DECIPHER
betapart
car
DescTools
iNEXT
aRtsy
vegan
ARTool
emmeans
lmerTest
patchwork
brms
ggvenn
maps
ggmap
sf
rgbif
rnaturalearth
paletteer
rentrez 

## ✅ Notes
- Processed files are generated within analysis scripts but may be reused across workflows  
- Scripts are designed to be run in sequence for reproducibility  
- Year-specific differences (2024 vs 2025) are handled explicitly in early processing

