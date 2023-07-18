#setwd("/Users/jeinhaus/Library/Mobile Documents/com~apple~CloudDocs/Documents/Stanford # #PostDoc/05 Drug Assay Preterm Birth/R analysis")
#funct_path = "/Users/jeinhaus/Library/Mobile Documents/com~apple~CloudDocs/Documents/Stanford PostDoc/05 #Drug Assay Preterm Birth/R analysis/Drug assay 15 drugs (Regated) - Statistics(1).csv"
#out_path = "/Users/jeinhaus/Library/Mobile Documents/com~apple~CloudDocs/Documents/Stanford PostDoc/05 Drug #Assay Preterm Birth/R analysis"

setwd("/Users/jonasamar/Stabl/Drug Study/Rscripts")
funct_path="../Drug Study/Drug assay csv/Drug assay ID1-4re-normalized&barcorded - Statistics - ID2P3 removed.csv"
out_path="../Drug Study/Drug assay csv"

## Libraries
library(tidyverse)
library(tidyselect)
library(reshape2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(paletteer)
library(tidyr)
library(janitor)
library(lme4)
library(lmerTest)

## Restructuring the Data
MyData = read.csv(funct_path)

# Removing redundant information (channel/reagent, uniquePopulationName/population) and useless columns (filename, parentPopulation)
MyData = subset(MyData, select = -c(filename, uniquePopulationName, parentPopulation, channel))

# Rewriting the Drug used on each Plate
#This vector contains the name of the drugs in the order corresponding to the plate number they were tested
Drugs <- c( 
  "Cefotaxime",
  "Lansoprazole",
  "Iopamidol",
  "Iohexol",
  "Benzylpenicillin",
  "Chlorthalidone",
  "Rifabutin",
  "Iodixanol",
  "Metformin",
  "Folic acid",
  "Clotrimazole",
  "Maprotiline",
  "Progesterone",
  "Pravastatin",
  "Methylpredonisolone"
)

# We associate the name of tje drug corresponding to its plate number in MyData
for (i in seq(1, 15)){
  MyData = MyData %>%
    mutate(Drug =  ifelse(Plate == i, Drugs[i], Drug))
}

# Removing the Plate column
MyData = subset(MyData, select = -c(Plate)) 

# Getting rid of the units for the doses and associating 0.5% to the null dose (it corresponds DMSO)
MyData = MyData %>%
  mutate(Dose =  ifelse(Dose %in% c("1ug", "1ng"), "1" , Dose))
MyData = MyData %>%
  mutate(Dose =  ifelse(Dose %in% c("10ug","10ng"), "10" , Dose))
MyData = MyData %>%
  mutate(Dose =  ifelse(Dose %in% c("100ug", "100ng"), "100" , Dose))
MyData = MyData %>%
  mutate(Dose =  ifelse(Dose %in% c("1000ug", "1000ng"), "1000" , Dose))
MyData = MyData %>%
  mutate(Dose =  ifelse(Dose %in% c("0.50%", "0.5%"), "0", Dose))

# Renaming the columns for convenience
MyData = MyData %>% rename(
  "dose" = "Dose",
  "drug" = "Drug",
  "stimulation" = "Stims"
)

## Preprocessing of the features

#Here use the asinh transform on the cytometry data and build a feature that is actually the difference between stimulated and un-stimulated samples .
final_data = c()
metadata = setdiff(colnames(MyData), c("stimulation", "median"))

for(dosepoint in unique(MyData$dose)){
  # filtering the data with a specific dose of drug
  MyData_dosepoint = MyData %>% filter(dose == dosepoint)
  # Getting the stims for this dose
  stims = unique(MyData_dosepoint$stimulation) 
  # asinh transformation
  MyData_dosepoint = MyData_dosepoint %>% mutate(feature = asinh(median/5))
  
  # Removing the median column and calculating Stimulated - Unstimulated
  MyData_dosepoint = MyData_dosepoint[, names(MyData_dosepoint) != "median"] %>%
    pivot_wider(names_from = stimulation, values_from = feature)
  
  for (stim in stims){
    if (stim!="Unstim"){
      MyData_dosepoint[stim] = MyData_dosepoint[stim] - MyData_dosepoint["Unstim"]
    }
  }
  
  MyData_fin = MyData_dosepoint %>% pivot_longer(-all_of(metadata),
                                                 names_to="stimulation",values_to="feature")
  
  # Adding the new calculated features to the final dataset
  final_data = rbind(final_data, MyData_fin)
}

write.csv(final_data, paste(out_path, "preprocessed.csv", sep = "/"))


## Data Analysis
#Here we want to quantify the behavior of the drugs into three categories : 'Inhibitor', 'Activation', 'FALSE'. Those categories are chosen based on the decrease or increase of the expression of each reagent when the drug dose increases.

### Approximating the behaviour of a drug on a specific feature
#Here we are creating a mixed linear models to quantify the effect of each drug on each feature.

Doseresponse = final_data %>%
  dplyr::mutate(dose = paste0("D", dose)) %>% # adds a D in front of the dose value
  dcast(ID + population + reagent + drug + stimulation ~ dose, value.var = "feature") %>% 
  dplyr::select(-D0) %>% # removing D0 (DMSO = no drug)
  mutate(feature=paste0(population, "_", reagent, "_", stimulation)) %>% # feature = population_reagent_stimulation
  pivot_longer(starts_with("D1")) %>%
  rename("dose" = "name") %>%
  select(c("ID", "feature", "drug", "dose", "value")) %>%
  na.omit() %>%
  mutate(dose = factor(dose),
         drug = factor(drug),
         feature = factor(feature),
         ID = factor(ID))

#Building a mixed linear regression on the entire dataset trying to identify the source of variances :
global.model <- lmer(value ~ dose + drug + feature + (1|drug:feature) + (1|ID), data=Doseresponse)