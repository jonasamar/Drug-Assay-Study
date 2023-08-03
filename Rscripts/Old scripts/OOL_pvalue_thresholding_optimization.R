# Looking for optimal pvalue thresholding for CYTOF + EGA

## Paths

#Copy Paste Here the path to your Drug-Assay-Study folder
my_path<-"/Users/jonasamar/Desktop/Drug-Assay-Study"

setwd(paste0(my_path,"/Rscripts"))
OOL_path=paste0(my_path,"/Onset of Labor csv")
drug_assay_path=paste0(my_path,"/Drug assay csv")

## Libraries

library(doSNOW)
library(parallel)
library(matrixStats)
library(zoo)
library(ggVennDiagram)
library(tidyverse)
library(readxl)
library(dplyr)
set.seed(2018)

## Importing and preprocessing data

# Preprocessing function which:
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
  cols <- colnames(prepro_data)
  # Imputing missing values with median of the colum
  for (col_name in cols[!cols %in% "ID"]){
    median <- median(prepro_data[[col_name]], na.rm = TRUE)
    prepro_data[[col_name]][is.na(prepro_data[[col_name]])] <- median
  }
  # Removing columns with standard deviation below std_threshold
  keep_cols <- which(colSds(as.matrix(prepro_data[,-which(names(prepro_data) == "ID")]))>std_threshold)
  prepro_data <- prepro_data[,names(keep_cols)]
  return(prepro_data)
}

#Importing and preprocessing the data for DOS between DOS_inf and DOS_sup

# DOS limits of the training data
DOS_inf=-120
DOS_sup=1
# all.simulated.data
load(paste0(OOL_path, "/immunome_noEGA_OOL_with_drugs.rda"))
# Outcomes
Yh <- read_csv(paste0(OOL_path, "/outcome_OOL.csv"), show_col_types = FALSE)
# CYTOF data
#CYTOF <- read_csv(paste0(OOL_path, "/immunome_noEGA_pen_OOL.csv"), show_col_types = FALSE)
CYTOF <- read_csv(paste0(OOL_path, "/immunome_EGA_pen_OOL.csv"), show_col_types = FALSE)
# Reorordering and filtering CYTOF data
data <- left_join(Yh, CYTOF, by="ID") %>% filter((DOS <= DOS_sup)&(DOS>=DOS_inf))
Yh <- data$DOS
CYTOF <- dplyr::select(data, -DOS)
# Preprocessing
CYTOF <- preprocessing(CYTOF, 0.2, 0.)
# Patients Id
Id <- as.factor(str_extract(data[["ID"]], "(?<=P)\\d+"))

## Model

# Helper function for predictive modeling and LOOCV
# Input arguments:
#   X: data matrix with patients, and omics measurments on rows and columns, respectively.
#   Y: response vector.
#   foldid: patients index that is a unique id for all subjects corresponding to the same patient.
#   i: the patient id that should be left out during the training process
#   parm: parameters of the Lasso 
# Output args,
#   ret: a list of predictions on the training and test data with the coefficients of Lasso 

#The complete analysis was performed using a similar strategy but using parallel processing to optimize lambda, followed by a stack generalization layer as described in the article.
xxx<-function(X, Y, foldid, i, parm, predict.simulated.data=FALSE)
{
  suppressMessages(library(randomForest, quietly = TRUE))
  suppressMessages(library(glmnet, quietly = TRUE))
  set.seed(2018+123*i)
  
  iInd=which(foldid==unique(foldid)[i])
  if(length(iInd)<2){
    iInd=c(iInd, iInd)
  }
  if(parm$scale=='ALL'){
    X=scale(X)
  }
  if(parm$scale=='Patient')
  {
    for(ap in seq(length(unique(foldid))))
    {
      sclidx= which(foldid==unique(foldid)[ap])
      if(length(sclidx)>1)
        X[sclidx]=scale(X[sclidx], scale=F)
    }
    X=scale(X) #creates NA values for 10 features
    rm_columns <- which(colSums(is.na(X)) > 0) # Getting the columns that are removed
    original.X=X # Keeping a copy of X without any columns removed
    X=X[, complete.cases(t(X))] #we remove those features
  }
  # Calculate Spearman's rank correlation and p-value
  results <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    res <- cor.test(X[, j], Y, method = "spearman", use = "pairwise.complete.obs", exact = FALSE)
    results[[j]] <- data.frame(colname = colnames(X)[j], pvalue = res$p.value)
  }
  
  # Combine results into a data frame
  pval_df <- bind_rows(results) %>%
    mutate(pvalue = as.numeric(pvalue))
  
  # Removing features with pvalue under a pval_threshold
  fiter.out.features <- pval_df[pval_df$pvalue >= parm$pval_threshold, "colname"]
  X <- X[,!(colnames(X) %in% fiter.out.features)]
  rm_columns <- c(setNames(match(fiter.out.features, colnames(original.X)), fiter.out.features), rm_columns)
  
  if (ncol(X) < 2){
    X=cbind(X,rep(1, nrow(X)))
  }
  
  XX=X[-iInd,]
  YY=Y[-iInd]
  XT=X[iInd,]
  fld=as.numeric(foldid[-iInd])
  
  # Transforming the values of fld so that they range from 1 to 52 (avoiding one warning from glmet)
  fld[fld>=i]=fld[fld>=i]-1
  
  ret = list()
  
  # Removed columns
  ret$rm_cols <- rm_columns
  
  # LASSO 1SE
  cvglm = cv.glmnet(XX, YY,  standardize=F, alpha=1, foldid = fld)
  ret$p1 = predict(cvglm, XT, s='lambda.1se')
  ret$ptrain1 = predict(cvglm, XX, s='lambda.1se')
  ret$coef = coef(cvglm, s='lambda.1se')[-1]
  ret$model = cvglm
  
  # Elastic Net
  EN_cv = cv.glmnet(x = XX, y = YY, foldid = fld, alpha = parm$a)
  ret$EN_p1 = predict(EN_cv, XT, s=EN_cv$lambda.min)
  ret$EN_ptrain1 = predict(EN_cv, XX, s=EN_cv$lambda.min)
  ret$EN_coef = coef(EN_cv, s=EN_cv$lambda.min)[-1]
  
  # Ridge Regression
  ridge_cv = cv.glmnet(x = XX, y = YY, foldid = fld, alpha = 0)
  ret$ridge_p1 = predict(ridge_cv, XT, s=ridge_cv$lambda.min)
  ret$ridge_ptrain1 = predict(ridge_cv, XX, s=ridge_cv$lambda.min)
  ret$ridge_coef = coef(ridge_cv, s=ridge_cv$lambda.min)[-1]
  
  # Adaptive LASSO
  alasso1_cv = cv.glmnet(x = XX, y = YY, foldid = fld, alpha = 1,
                         penalty.factor = 1 / abs(ret$ridge_coef),
                         keep = TRUE)
  ret$alasso_p1 = predict(alasso1_cv, XT, s=alasso1_cv$lambda.min)
  ret$alasso_ptrain1 = predict(alasso1_cv, XX, s=alasso1_cv$lambda.min)
  ret$alasso_coef = coef(alasso1_cv, s=alasso1_cv$lambda.min)[-1]
  
  # Random Forest
  rd_forest_cv = randomForest(x = XX, y = YY, importance=TRUE)
  ret$rd_forest_p1 = predict(rd_forest_cv, XT)
  ret$rd_forest_ptrain1 = predict(rd_forest_cv, XX)
  ret$rd_forest_importance = rd_forest_cv$importance
  
  if (predict.simulated.data){
    # Predictions of Adaptive LASSO on simulated data
    for (target_drug in unique(all.simulated.data$drug)){
      # Applying CYTOF preprocessing to the simulated data
      prepro_data <- all.simulated.data %>%
        # Filtering the samples of the patients in the test set
        filter(drug==target_drug,
               grepl(paste("P",levels(foldid)[i],"_",sep=""),ID)) %>%
        # Removing the columns which were removed in CYTOF
        dplyr::select(colnames(CYTOF), -names(rm_columns)) %>%
        # Imputing the missing values with the medians of the corresponding CYTOF columns
        mutate(across(everything(), ~ ifelse(is.na(.), median(CYTOF[[cur_column()]], na.rm = TRUE), .))) %>%
        # Converting to a matrix
        column_to_rownames(var = "ID") %>%
        as.matrix()
      
      # Predictions
      if (parm$predict.model == "Random.Forest"){
        ret[[target_drug]] <- predict(rd_forest_cv, prepro_data)
      }
      if (parm$predict.model == "ElasticNet"){
        ret[[target_drug]] <- predict(EN_cv, prepro_data, s=EN_cv$lambda.min)
      }
      if (parm$predict.model == "Lasso.1se"){
        ret[[target_drug]] <- predict(cvglm, prepro_data, s='lambda.1se')
      }
    }
  }
  
  return(ret)
}

#Prediction of the model
# Function which takes the outcome of xxx and aggregate the predictions across all the models
aggregation.of.predictions <- function(prdC){
  list_name_pred.key <- list("Lasso.1se"="p1", 
                             "Ridge.Regression"="ridge_p1",
                             "ElasticNet"="EN_p1",
                             "Adaptive.Lasso"="alasso_p1",
                             "Random.Forest"="rd_forest_p1")
  list_preds <- list()
  for (model.name in names(list_name_pred.key)){
    pred.key=list_name_pred.key[[model.name]]
    ccc=vector()
    for(i in seq(npt))
    {
      iInd=which(Id==unique(Id)[i])
      if(length(iInd)>1)
      {
        ccc[iInd] = prdC[[i]][[pred.key]]
      }
      else
      {
        ccc[iInd] = prdC[[i]][[pred.key]][1]
      }
    }
    list_preds[[model.name]]=ccc
  }
  return(list_preds)
}

## Optimization the p-value threshold

#Here we optimize the pvalue_threshold to maximize the spearmanr score of our final model

# Range of thresholds tested
pval_thresholds=seq(1e-14, 1e-13, length.out=51)

# xxx parameters
parm=list()
parm$scale='Patient'
parm$a=0.5
npt=length(unique(Id))

# Initializing the list of the spearamnr scores
listr<-c()

# Filling listr
for (pval_threshold in pval_thresholds){
  parm$pval_threshold=pval_threshold
  print(pval_threshold)
  # Opening backend-configuration
  cl <- makeSOCKcluster(detectCores())
  registerDoSNOW(cl)
  # Getting the predictions
  prdC=foreach(i=seq(npt),.packages = c("dplyr", "tibble")) %dopar% xxx(data.matrix(CYTOF), Yh, Id, i, parm)
  # Closing backend configuration
  stopCluster(cl)
  # Getting the list of the best spearmanr score among the 5 models
  list_preds <- aggregation.of.predictions(prdC)
  # Extracting the maximum spearmanr score
  spearmr <- c()
  for (model.name in names(list_preds)){
    y_pred <- list_preds[[model.name]]
    spearmr <- c(spearmr, cor(Yh, y_pred, method = "spearman"))
  }
  # Adding the best spearmanr score to the list
  listr <- rbind(listr, c(pval_threshold, length(prdC[[1]]$rm_cols), spearmr, max(spearmr)))
}
# Renaming the columns
colnames(listr) <- c("pval_threshold", "n_rm_cols", names(list_preds), "best_r_score")
# Saving the dataframe
save(listr, file=(paste0(OOL_path, "/Pvalue threshold optimization with EGA.rda")))