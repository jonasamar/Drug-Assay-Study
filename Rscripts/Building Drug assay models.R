## Paths

#Copy Paste Here the path to your Drug-Assay-Study folder
my_path<-"/Users/jonasamar/Desktop/Drug-Assay-Study"

setwd(paste0(my_path,"/Rscripts"))
out_path=paste0(my_path,"/Drug assay csv")

## Libraries

library(tidyverse)
library(tidyselect)
library(reshape2)
library(lme4)
#library(lmerTest)
library(ggplot2)
#library(r2glmm)
library(sjstats)
library(doSNOW)
library(parallel)
library(reticulate)
library(aomisc)
library(dplyr)
library(dplyr)

## Building model function

# Function which attempts to build sigmoid models on the data
build.sigmoid.or.linear.model <- function(data, feature.name, drug.name, id="median"){
  # Building a linear or sigmoid model only when the response is not constant (we get an error message otherwise...)
  if (length(unique(data$value)) > 1){
    # Generate a sequence of intermediate doses
    intermediate_dose_0 <- lapply(1:(length(data$dose) - 1), function(i) (data$dose[i + 1] + data$dose[i])/2)
    intermediate_dose_1 <- lapply(1:(length(data$dose) - 1), function(i) data$dose[i]+(data$dose[i + 1] - data$dose[i])/3)
    intermediate_dose_2 <- lapply(1:(length(data$dose) - 1), function(i) data$dose[i]+2*(data$dose[i + 1] - data$dose[i])/3)
    intermediate_dose_3 <- lapply(1:(length(data$dose) - 1), function(i) data$dose[i]+(data$dose[i + 1] - data$dose[i])/4)
    intermediate_dose_4 <- lapply(1:(length(data$dose) - 1), function(i) data$dose[i]+3*(data$dose[i + 1] - data$dose[i])/4)
    intermediate_dose <- c(unlist(intermediate_dose_0), unlist(intermediate_dose_1), unlist(intermediate_dose_2), unlist(intermediate_dose_3), unlist(intermediate_dose_4), data$dose)
    # Interpolate the intermediate values at the intermediate doses
    interpolated_values_0 <- lapply(1:(length(data$value) - 1), function(i) (data$value[i + 1] + data$value[i])/2)
    interpolated_values_1 <- lapply(1:(length(data$value) - 1), function(i) data$value[i]+2*(data$value[i + 1] - data$value[i])/100)
    interpolated_values_2 <- lapply(1:(length(data$value) - 1), function(i) data$value[i]+98*(data$value[i + 1] - data$value[i])/100)
    interpolated_values_3 <- lapply(1:(length(data$value) - 1), function(i) data$value[i]+(data$value[i + 1] - data$value[i])/100)
    interpolated_values_4 <- lapply(1:(length(data$value) - 1), function(i) data$value[i]+99*(data$value[i + 1] - data$value[i])/100)
    interpolated_values <- c(unlist(interpolated_values_0), unlist(interpolated_values_1), unlist(interpolated_values_2),unlist(interpolated_values_3), unlist(interpolated_values_4), data$value)
    # Create a data frame with the interpolated data
    interpolated_data <- data.frame(dose = intermediate_dose, value = interpolated_values)
    # Sort the interpolated data by dose
    interpolated_data <- interpolated_data[order(interpolated_data$dose), ]
    
    # Sigmoid function we want to try to fit
    sigmoid_model <- function(x, b, c, d, e) {
      c + (d - c) / (1 + exp(-b * (x - e)))
    }
    model <- NULL
    class <- NULL
    # Attempt to fit the sigmoid model using nls() and fall back to a linear model
    try_result <- tryCatch({
      #Initilize parameters
      gaps <- abs(diff(data$value))
      jump <- which.max(gaps)
      start_e=mean(c(data$dose[jump], data$dose[jump+1]))
      start_b=sign(data$value[jump+1]-data$value[jump])
      startp <- list(b=start_b, c=min(data$value), d=max(data$value), e=start_e)
      # Model fitting
      model <-nls(value ~ sigmoid_model(dose, b, c, d, e), data = interpolated_data, start=startp, algorithm = "port")
      class <- "sigmoid"
    }, error = function(e) {
      model<<-suppressMessages(lm(value ~ dose, data = data))
      class <<- "linear"
    })
    
    # Building scores
    if (class=="sigmoid") {
      # Sigmoid model fitting was successful
      coefs <- coef(model)
      preds <- sigmoid_model(data$dose, coefs[["b"]], coefs[["c"]], coefs[["d"]], coefs[["e"]])
      residuals <- preds - data$value
      rmse <- sqrt(mean(residuals^2))
      aic <- 8 - 2*log(1/sqrt(2*pi*sum(residuals^2)/nrow(interpolated_data))*exp(-sum(residuals^2)/(2*nrow(interpolated_data))))
      # Trying to get a pvalue when available
      pval<-NULL
      try_pval <- tryCatch({
        pval<-summary(model)$coefficients[, "Pr(>|t|)"][["b"]]
      }, error = function(e) {
        pval<<-NA
      })
    }
    if (class=="linear"){
      # Sigmoid model fitting was NOT successful
      residuals <- residuals(model)
      rmse <- sqrt(mean(residuals^2))
      aic <- AIC(model)
      # Trying to get a pvalue when available
      pval<-NULL
      try_pval <- tryCatch({
        pval <- summary(model)$coefficients[, "Pr(>|t|)"][["dose"]]
      }, error = function(e) {
        pval<<-NA
      })
    }
  }else{
    # When there is only one value, the model is constant
    model <- list(unique(data$value))
    aic <- NA
    pval <- NA
    rmse <- NA
    class <- "constant"
  }
  # Adding new linear model to models
  return(data.frame(feature = feature.name,
                    drug = drug.name,
                    ID=id,
                    model = I(list(model)),
                    aic=aic,
                    rmse=rmse,
                    pval=pval,
                    class = class,
                    stringsAsFactors = FALSE))
}

## Individual models

#Importing penalized preprocessed data with the different scales.

# AllPenDoseresponse
load(paste0(out_path, "/all pen dose response with scales.rda"))

#Function building a the mixed models for a given drug across the 720 features.
build.individual.model <- function(drug.name, dataset){
  # Create an empty data frame to store the results
  models <- data.frame(stringsAsFactors = FALSE)
  
  for (feature.name in unique(dataset$feature)){
    # Subset of the dataset 
    Subset <- dataset %>% 
      filter(feature == feature.name) %>%
      filter(drug == drug.name)  %>%
      dplyr::select(c("ID", "dose", "value"))
    for (id in unique(Subset$ID)){
      # Subset at the individual level
      subset <- Subset %>%
        filter(ID==id) %>%
        dplyr::select(c("ID", "dose", "value"))
      new_model <- build.sigmoid.or.linear.model(subset,feature.name, drug.name,id)
      # Adding new linear model to models
      models <- rbind(models,new_model)
    }
  }
  return(models)
}

#log10.individual.model : dataframe containing all the informations about the mixed linear models that are built on each feature on a log10 scale of the concentrations.
# Setting the dose values to their log10 transform
Doseresponse <- AllPenDoseresponse %>%
  mutate(dose = log10_dose) %>%
  filter(is.finite(dose))

# Set up parallel backend using doParallel
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

# Create individual linear models on log10 scale
log10.individual.model = foreach(drug.name = unique(Doseresponse$drug), 
                                 .packages = c("tidyverse","lme4", "lmerTest","dplyr"),
                                 .combine = rbind) %dopar% build.individual.model(drug.name, Doseresponse)

# Close parallel backend
stopCluster(cl)

# Saving mixed models
save(log10.individual.model, 
     file=paste0(out_path, "/log10 individual models.rda"))

#log10_1.individual.model : dataframe containing all the informations about the mixed linear models that are built on each feature on a log10(x+1) scale of the concentrations.
# Setting the dose values to their log10(x+1) transform
Doseresponse <- AllPenDoseresponse %>%
  mutate(dose = log10_1_dose) %>%
  filter(is.finite(dose))

# Set up parallel backend using doParallel
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

# Create individual linear models on log10 scale
log10_1.individual.model = foreach(drug.name = unique(Doseresponse$drug), 
                                   .packages = c("tidyverse","lme4", "lmerTest", "dplyr"),
                                   .combine = rbind) %dopar% build.individual.model(drug.name, Doseresponse)

# Close parallel backend
stopCluster(cl)

# Saving mixed models
save(log10_1.individual.model, 
     file=paste0(out_path, "/log10_1 individual models.rda"))

#pseudo_log10.individual.model : dataframe containing all the informations about the mixed linear models that are built on each feature on a pseudo_log10 (=asinh(x/2)/log(10)) scale of the concentrations.
# Setting the dose values to their pseudo log transform
Doseresponse <- AllPenDoseresponse %>%
  mutate(dose = pseudo_log_dose) %>%
  filter(is.finite(dose))

# Set up parallel backend using doParallel
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

# Create linear mixed models on pseudo log scale
pseudo_log10.individual.model = foreach(drug.name = unique(Doseresponse$drug), 
                                        .packages = c("tidyverse","lme4", "lmerTest", "dplyr"),
                                        .combine = rbind) %dopar% build.individual.model(drug.name, Doseresponse)

# Close parallel backend
stopCluster(cl)

# Saving mixed models
save(pseudo_log10.individual.model, 
     file=paste0(out_path, "/pseudo_log10 individual models.rda"))


## Linear models on medians

#Importing penalized preprocessed data with the different scales.

# MedianDoseresponse
load(paste0(out_path, "/median pen dose response with scales.rda"))

#Function building a the linear models for a given drug across the 720 features.
build.median.models <- function(drug.name, dataset){
  # Create an empty data frame to store the results
  models <- data.frame(stringsAsFactors = FALSE)
  
  # We look at each drug for a given feature
  for (feature.name in unique(dataset$feature)){
    # subset of the dataset 
    subset <- dataset %>% 
      filter(feature == feature.name) %>%
      filter(drug == drug.name)  %>%
      mutate(value=median_value) %>%
      dplyr::select(c("dose", "value")) %>%
      distinct() %>%
      suppressMessages()
    
    new_model <- build.sigmoid.or.linear.model(subset,feature.name, drug.name)
    # Adding new linear model to models
    models <- rbind(models,new_model)
  }
  return(models)
}

#log10.median.model : dataframe containing all the information about the linear models that are built on each feature on a log10 scale of the concentrations.
# Setting the dose values to their log10 transform
Doseresponse <- MedianDoseresponse %>%
  mutate(dose = log10_dose) %>%
  filter(is.finite(dose))

# Set up parallel backend using doParallel
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

# Create linear mixed models on log10 scale
log10.median.model = foreach(drug.name = unique(Doseresponse$drug), 
                             .packages = c("tidyverse","lme4", "lmerTest", "dplyr"),
                             .combine = rbind) %dopar% build.median.models(drug.name, Doseresponse)

# Close parallel backend
stopCluster(cl)

# Saving mixed models
save(log10.median.model, 
     file=paste0(out_path, "/log10 median models.rda"))

#log10_1.median.model : dataframe containing all the information about the linear models that are built on each feature on a log10(x+1) scale of the concentrations.
# Setting the dose values to their log10 transform
Doseresponse <- MedianDoseresponse %>%
  mutate(dose = log10_1_dose) %>%
  filter(is.finite(dose))

# Set up parallel backend using doParallel
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

# Create linear mixed models on log10 scale
log10_1.median.model = foreach(drug.name = unique(Doseresponse$drug),
                               .packages = c("tidyverse","lme4", "lmerTest", "dplyr"),
                               .combine = rbind) %dopar% build.median.models(drug.name, Doseresponse)

# Close parallel backend
stopCluster(cl)

# Saving mixed models
save(log10_1.median.model, 
     file=paste0(out_path, "/log10_1 median models.rda"))

#pseudo_log10.median.model : dataframe containing all the information about the linear models that are built on each feature on a pseudo_log10 scale of the concentrations.
# Setting the dose values to their log10 transform
Doseresponse <- MedianDoseresponse %>%
  mutate(dose = pseudo_log_dose) %>%
  filter(is.finite(dose))

# Set up parallel backend using doParallel
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

# Create linear mixed models on log10 scale
pseudo_log10.median.model = foreach(drug.name = unique(Doseresponse$drug),
                                    .packages = c("tidyverse","lme4", "lmerTest", "dplyr"),
                                    .combine = rbind) %dopar% build.median.models(drug.name, Doseresponse)

# Close parallel backend
stopCluster(cl)

# Saving mixed models
save(pseudo_log10.median.model, 
     file=paste0(out_path, "/pseudo_log10 median models.rda"))