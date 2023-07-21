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
library(gridExtra)
library(grid)
set.seed(2018)

## Importing and preprocessing data

# Onset of Labor data
load(paste0(OOL_path, "/OOL paper data/Data.rda"))
##All features not consistent with the prior immunological knowledge were excluded from the modeling process.
CPen = cbind(CPen, matrix(1, nrow = nrow(CPen) , ncol = 41))
CYTOF = CYTOF[,which(CPen[1,]!=0)]
CSds = which(colSds(CYTOF)!=0)
Proteomics = Proteomics
Metabolomics = Metabolomics
# Adding EGA
CYTOF = cbind(EGA, CYTOF)
Proteomics = cbind(EGA, Proteomics)
Metabolomics = cbind(EGA, Metabolomics)
# Getting the preterm rows
preterm_rows <- which(grepl(paste(c("017", "008", "003", "005", "027"), collapse = "|"), rownames(CYTOF)))

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
pval_thresholds=c(seq(1e-13, 1e-10, length.out=51), seq(1e-10, 5e-2, length.out=51))

# xxx parameters
parm=list()
parm$scale='Patient'
parm$a=0.5
npt=length(unique(Id))

# Initializing the list of the spearamnr scores
listr <- list("CYTOF"=list(), "Proteomics"=list(), "Metabolomics"=list())

# Filling listr
for (omic in names(listr)){
  listr_omic = list()
  for (pval_threshold in pval_thresholds){
    parm$pval_threshold=pval_threshold
    print(pval_threshold)
    # Opening backend-configuration
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    # Getting the predictions
    prdC=foreach(i=seq(npt),.packages = c("dplyr", "tibble")) %dopar% xxx(data.matrix(CYTOF), DOS, Id, i, parm)
    # Closing backend configuration
    stopCluster(cl)
    # Getting the list of the best spearmanr score among the 5 models
    list_preds <- aggregation.of.predictions(prdC)
    # Extracting the maximum spearmanr score
    spearmr <- c()
    for (model.name in names(list_preds))
      y_pred <- list_preds[[model.name]]
      spearmr <- c(spearmr, cor(DOS, y_pred, method = "spearman"))
    }
    # Adding the best spearmanr score to the list
    listr_omic <- rbind(listr_omic, c(pval_threshold, length(prdC[[1]]$rm_cols), spearmr, max(spearmr)))
  }
  # Renaming the columns
  colnames(listr_omic) <- c("pval_threshold", "n_rm_cols", names(list_preds), "best_r_score")
  listr[[omic]]=listr_omic
}

# Saving the dataframe
save(listr, file=(paste0(OOL_path, "/Pvalue threshold optimization for OOL data with EGA.rda")))

## Model fitting

#CYTOF
# xxx parameter
parm=list()
parm$scale='Patient'
parm$a=0.5
npt=length(unique(Id))
parm$pval_threshold=listr$CYTOF$pval_threshold[which.max(listr$CYTOF$speamanr)]
# Opening backend-configuration
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)
# Getting models and predictions
prdC=foreach(i=seq(npt),.packages = c("dplyr", "tibble")) %dopar% xxx(data.matrix(CYTOF), DOS, Id, i, parm, FALSE) 
list_preds_C<-aggregation.of.predictions(prdC)
# Closing backend configuration
stopCluster(cl)

#Proteomics
# xxx parameter
parm=list()
parm$scale='Patient'
parm$a=0.5
npt=length(unique(Id))
parm$pval_threshold=listr$Proteomics$pval_threshold[which.max(listr$Proteomics$speamanr)]

# Opening backend-configuration
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)
# Getting models and predictions
prdP=foreach(i=seq(npt),.packages = c("dplyr", "tibble")) %dopar% xxx(data.matrix(Proteomics), DOS, Id, i, parm, FALSE) 
list_preds_P<-aggregation.of.predictions(prdP)
# Closing backend configuration
stopCluster(cl)


#Metabolomics
# xxx parameter
parm=list()
parm$scale='Patient'
parm$a=0.5
npt=length(unique(Id))
parm$pval_threshold=listr$Metabolomics$pval_threshold[which.max(listr$Metabolomics$speamanr)]

# Opening backend-configuration
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)
# Getting models and predictions
prdM=foreach(i=seq(npt),.packages = c("dplyr", "tibble")) %dopar% xxx(data.matrix(Metabolomics), DOS, Id, i, parm, FALSE) 
list_preds_M<-aggregation.of.predictions(prdM)
# Closing backend configuration
stopCluster(cl)

## Evaluation of results

#Performances of the models
require(Metrics)

# Plot function
plot.model.results <- function(y_pred, y_true, title, x_name="DOS"){
  myPv = cor.test(y_true, y_pred, method = 'spearman', exact = FALSE)$p.value
  myerr = sqrt(mean((y_true-y_pred)^2))
  myr2 = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  mycorr = cor(y_true, y_pred, method = "spearman")
  
  if (x_name=="DOS"){dash_coefs=c(0.,1.)}
  if (x_name=="EGA"){dash_coefs=c(-259, 6.5)}
  
  data <- data.frame(x = y_true, y = y_pred)
  p <- ggplot(data, aes(x = x, y = y)) + 
    geom_point() +
    geom_abline(intercept = dash_coefs[1], slope = dash_coefs[2], linetype = "dashed") +
    labs(x = x_name, y = "Prediction") +
    ggtitle(paste("Model : ", title,
                  "\nSpearmanr : ", round(mycorr, digits = 5), 
                  "\nRMSE : ",round(myerr, digits = 5), 
                  "\np-value : ",round(myPv, digits = 30), 
                  "\nR^2 : ",round(myr2, digits = 5), sep=""))
  print(p)
}

# Plotting the results of the models for CYTOF omic
plots <- list()
for (model.name in names(list_preds_C)){
  # with DOS
  plots[[paste0(model.name,1)]] <- plot.model.results(list_preds_C[[model.name]][preterm_rows], DOS[preterm_rows], paste0(model.name, " on preterm"), "DOS")
  plots[[paste0(model.name,2)]] <- plot.model.results(list_preds_C[[model.name]][-preterm_rows], DOS[-preterm_rows], paste0(model.name, " on term"), "DOS")
  plots[[paste0(model.name,3)]] <- plot.model.results(list_preds_C[[model.name]], DOS, paste0(model.name, " on all"), "DOS")
  # with EGA
  plots[[paste0(model.name,4)]] <- plot.model.results(list_preds_C[[model.name]][preterm_rows], EGA[preterm_rows], paste0(model.name, " on preterm"), "EGA")
  plots[[paste0(model.name,5)]] <- plot.model.results(list_preds_C[[model.name]][-preterm_rows], EGA[-preterm_rows], paste0(model.name, " on term"), "EGA")
  plots[[paste0(model.name,6)]] <- plot.model.results(list_preds_C[[model.name]], EGA, paste0(model.name, " on all"), "EGA")
}
pdf(file=paste0(OOL_path, "/plots/OOL paper CYTOF with EGA best models.pdf"), height=8, width=8)
for (i in 0:floor(length(plots)/6)){
  grid.arrange(grobs = plots[(6*i+1):(6*i+6)], nrow = 2, ncol = 3)
  grid.newpage()
}
dev.off()

# Plotting the results of the models for Proteomics omic
plots <- list()
for (model.name in names(list_preds_P)){
  # with DOS
  plots[[paste0(model.name,1)]] <- plot.model.results(list_preds_P[[model.name]][preterm_rows], DOS[preterm_rows], paste0(model.name, " on preterm"), "DOS")
  plots[[paste0(model.name,2)]] <- plot.model.results(list_preds_P[[model.name]][-preterm_rows], DOS[-preterm_rows], paste0(model.name, " on term"), "DOS")
  plots[[paste0(model.name,3)]] <- plot.model.results(list_preds_P[[model.name]], DOS, paste0(model.name, " on all"), "DOS")
  # with EGA
  plots[[paste0(model.name,4)]] <- plot.model.results(list_preds_P[[model.name]][preterm_rows], EGA[preterm_rows], paste0(model.name, " on preterm"), "EGA")
  plots[[paste0(model.name,5)]] <- plot.model.results(list_preds_P[[model.name]][-preterm_rows], EGA[-preterm_rows], paste0(model.name, " on term"), "EGA")
  plots[[paste0(model.name,6)]] <- plot.model.results(list_preds_P[[model.name]], EGA, paste0(model.name, " on all"), "EGA")
}
pdf(file=paste0(OOL_path, "/plots/OOL paper Proteomics with EGA best models.pdf"))
for (i in 0:floor(length(plots)/6)){
  grid.arrange(grobs = plots[(6*i+1):(6*i+6)], nrow = 2, ncol = 3)
  grid.newpage()
}
dev.off()

# Plotting the results of the models for Metabolomics omic
plots <- list()
for (model.name in names(list_preds_M)){
  # with DOS
  plots[[paste0(model.name,1)]] <- plot.model.results(list_preds_M[[model.name]][preterm_rows], DOS[preterm_rows], paste0(model.name, " on preterm"), "DOS")
  plots[[paste0(model.name,2)]] <- plot.model.results(list_preds_M[[model.name]][-preterm_rows], DOS[-preterm_rows], paste0(model.name, " on term"), "DOS")
  plots[[paste0(model.name,3)]] <- plot.model.results(list_preds_M[[model.name]], DOS, paste0(model.name, " on all"), "DOS")
  # with EGA
  plots[[paste0(model.name,4)]] <- plot.model.results(list_preds_M[[model.name]][preterm_rows], EGA[preterm_rows], paste0(model.name, " on preterm"), "EGA")
  plots[[paste0(model.name,5)]] <- plot.model.results(list_preds_M[[model.name]][-preterm_rows], EGA[-preterm_rows], paste0(model.name, " on term"), "EGA")
  plots[[paste0(model.name,6)]] <- plot.model.results(list_preds_M[[model.name]], EGA, paste0(model.name, " on all"), "EGA")
}
pdf(file=paste0(OOL_path, "/plots/OOL paper Metabolomics with EGA best models.pdf"))
for (i in 0:floor(length(plots)/6)){
  grid.arrange(grobs = plots[(6*i+1):(6*i+6)], nrow = 2, ncol = 3)
  grid.newpage()
}
dev.off()