## Paths
setwd("~/OOL_model_opti")
OOL_path="~/OOL_model_opti/Onset of Labor csv"
drug_assay_path="~/OOL_model_opti/Drug assay csv"

## Libraries

library(doSNOW)
library(parallel)
library(matrixStats)
library(zoo)
library(ggVennDiagram)
library(tidyverse)
library(readxl)
set.seed(2018)

## Importing and preprocessing data

#Preprocessing function which:
#  1 - Remove features with proportion of NA above a certain threshold
#  2 - Replace remaining NA values with the median of the feature
#  3 - Remove features with standard deviation below a certain threshold
preprocessing <- function(dataset, NA_threshold, std_threshold){
  prepro_data <- dataset
  cols <- colnames(dataset)
  # Removing columns with proportion of missing values above NA_threshold
  rm_cols <- c()
  for (col_name in cols[!cols %in% "ID"]){
    na_count <- sum(is.na(prepro_data[[col_name]]))
    if (na_count/length(prepro_data[[col_name]]) > NA_threshold){
      rm_cols[length(rm_cols)+1] <- col_name
    }
  }
  prepro_data <- prepro_data[,!colnames(prepro_data) %in% rm_cols]
  # Imputing missing values with median of the colum
  for (col_name in cols[!cols %in% "ID"]){
    median <- median(prepro_data[[col_name]], na.rm = TRUE)
    prepro_data[[col_name]][is.na(prepro_data[[col_name]])] <- median
  }
  # Removing columns with standard deviation below std_threshold
  keep_cols <- which(colSds(as.matrix(prepro_data[,-which(names(prepro_data) == "ID")]))>std_threshold)
  prepro_data <- prepro_data[,keep_cols]
  return(prepro_data)
}

#Importing and preprocessing the data

# all.simulated.data
load(paste0(OOL_path, "/immunome_noEGA_OOL_with_drugs.rda"))
# Outcomes
Yh <- read_csv(paste0(OOL_path, "/outcome_OOL.csv"), show_col_types = FALSE)
# CYTOF data
CYTOF <- read_csv(paste0(OOL_path, "/immunome_noEGA_pen_OOL.csv"), show_col_types = FALSE)
# Reorordering Yh
Yh <- Yh[match(CYTOF$ID, Yh$ID), ][["DOS"]]
# Preprocessing
CYTOF <- preprocessing(CYTOF, 0.2, 0.)
# Patients Id
Id <- as.factor(str_extract(CYTOF[["ID"]], "(?<=P)\\d+"))
