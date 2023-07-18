library(doSNOW)
library(parallel)
library(matrixStats)
library(readxl)
set.seed(2018)
##The data file contains the following variables:
##CPen: Cell signaling-based penalization matrix.
##DOS: Days to onset of labor         
##EGA: Estimated Gestational Age         
##Timepoint: trimesters of pregnancy from G1 to G3 that represent the first three trimesters from the beginning to onset of labor.
##Id: patients index.
##CYTOF, Metabolomics, Proteomics: constitute the 3 omics datasets.

paper_path <- '~/Desktop/OOL_Code/Data.rda'
my_path <- '/Users/jonasamar/Stabl/Onset of Labor/Data/Data.rda'
  
load(my_path)
Yh = DOS

##All features not consistent with the prior immunological knowledge were excluded from the modeling process.
CPen = cbind(CPen, matrix(1, nrow = nrow(CPen) , ncol = 41))
CYTOF = CYTOF[,which(CPen[1,]!=0)]
CSds = which(colSds(CYTOF)!=0)
Proteomics = Proteomics
Metabolomics = Metabolomics
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

##Helper function for predictive modeling and LOOCV
## Input arguments:
##X: data matrix with patients, and omics measurments on rows and columns, respectively.
##Y: response vector.
##foldid: patients index that is a unique id for all subjects corresponding to the same patient.
##i: the patient id that should be left out during the training process
##parm: parameters of the Lasso 
##Output args,
##ret: a list of predictions on the training and test data with the coefficients of Lasso 
##The complete analysis was performed using a similar strategy but using parallel processing to optimize lambda, followed by a stack generalization layer as described in the article.
xxx<-function(X, Y, foldid, i, parm)
{
  suppressMessages(library(randomForest, quietly = TRUE))
  suppressMessages(library(glmnet, quietly = TRUE))
  set.seed(2018+123*i)
  
  iInd=which(foldid==unique(foldid)[i])
  if(length(iInd)<2)
    iInd=c(iInd, iInd)
  if(parm$scale=='ALL')
    X=scale(X)
  if(parm$scale=='Patient')
  {
    for(ap in seq(length(unique(foldid))))
    {
      sclidx= which(foldid==unique(foldid)[ap])
      if(length(sclidx)>1)
        X[sclidx]=scale(X[sclidx], scale=F)
    }
    X=scale(X)
  }
  
  XX=X[-iInd,]
  YY=Y[-iInd]
  XT=X[iInd,]
  fld=as.numeric(foldid[-iInd])
  
  ret = list()
  cvglm = cv.glmnet(XX, YY,  standardize=F, alpha=parm$a, foldid = fld)
  ret$p1 = predict(cvglm, XT, s='lambda.1se')
  ret$ptrain1 = predict(cvglm, XX, s='lambda.1se')
  ret$coef = coef(cvglm, s='lambda.1se')[-1]
  
  return(ret)
}
##Model fitting.
parm=list()
parm$scale='Patient'
parm$a=1
npt=length(unique(Id))
prdP=foreach(i=seq(npt)) %dopar% xxx(Proteomics, Yh, Id, i, parm)
prdM=foreach(i=seq(npt)) %dopar% xxx(Metabolomics, Yh, Id, i, parm)
prdC=foreach(i=seq(npt)) %dopar% xxx(CYTOF[,CSds], Yh, Id, i, parm)

##Evaluation of results.
myPvs=vector()

##proteomics
ppp=vector()
for(i in seq(npt))
{
  iInd=which(Id==unique(Id)[i])
  if(length(iInd)>1)
  {
    ppp[iInd] = prdP[[i]]$p1
  }
  else
  {
    ppp[iInd] = prdP[[i]]$p1[1]
  }
}
myPvs[3]=-log10(cor.test(Yh, ppp, method = 'spearman', exact = FALSE)$p.value)

##metabolomics
mmm=vector()
for(i in seq(npt))
{
  iInd=which(Id==unique(Id)[i])
  if(length(iInd)>1)
  {
    mmm[iInd] = prdM[[i]]$p1
  }
  else
  {
    mmm[iInd] = prdM[[i]]$p1[1]
  }
}
myPvs[2]=-log10(cor.test(Yh, mmm, method = 'spearman', exact = FALSE)$p.value)

##Immune-System (CyTOF)
ccc=vector()
for(i in seq(npt))
{
  iInd=which(Id==unique(Id)[i])
  if(length(iInd)>1)
  {
    ccc[iInd] = prdC[[i]]$p1
  }
  else
  {
    ccc[iInd] = prdC[[i]]$p1[1]
  }
}
myPvs[1]=-log10(cor.test(Yh, ccc, method = 'spearman', exact = FALSE)$p.value)

##Stacked generalization
SG=cbind(mmm, ppp, ccc)
prdSG=foreach(i=seq(npt)) %dopar% xxx(SG, Yh, Id, i, parm)

sss=vector()
for(i in seq(npt))
{
  iInd=which(Id==unique(Id)[i])
  if(length(iInd)>1)
  {
    sss[iInd] = prdSG[[i]]$p1
  }
  else
  {
    sss[iInd] = prdSG[[i]]$p1[1]
  }
}

myPvs[4]=-log10(cor.test(Yh, sss, method = 'spearman', exact = FALSE)$p.value)

##Model contribution plot in terms of -log10(p-value)
mydatasets=c('ImmuneSystem','Metabolomics', 'PlasmaSomalogic', 'SG')
mycols=c('#AB47BC', '#3F51B5','#F57C00', 'black')
par(mar=c(10,5,5,0))
barplot(myPvs, las=2, ylab='Model Contribution', col= mycols, names.arg = mydatasets, beside=T, ylim = c(0, 40))
abline(h=-log10(0.05), col="red", lty=2, lwd=3)

##Model contribution plot in terms of "Root Mean Square Error"
require(Metrics)
myerr = c(rmse(ccc, Yh), rmse(mmm, Yh), rmse(ppp, Yh), rmse(sss, Yh))
par(mar=c(10,5,5,0))
barplot(myerr, las=2, ylab='RMSE', col= mycols, names.arg = mydatasets, beside=T)

##Model contribution plot in terms of "R^2"
myr2 = c(1 - sum((Yh - ccc)^2) / sum((Yh - mean(Yh))^2),
         1 - sum((Yh - mmm)^2) / sum((Yh - mean(Yh))^2),
         1 - sum((Yh - ppp)^2) / sum((Yh - mean(Yh))^2),
         1 - sum((Yh - sss)^2) / sum((Yh - mean(Yh))^2))
par(mar=c(10, 5, 5, 0))
barplot(myr2, las=2, ylab='R-squared', col=mycols, names.arg=mydatasets, beside=T)

##Model contribution plot in terms of "Spearman R correlation"
mycorr = c(cor(Yh, ccc, method = "spearman"),
           cor(Yh, mmm, method = "spearman"),
           cor(Yh, ppp, method = "spearman"),
           cor(Yh, sss, method = "spearman"))

par(mar = c(10, 5, 5, 0))
barplot(mycorr, las = 2, ylab = "Spearman Correlation", col = mycols, names.arg = mydatasets, beside = TRUE)
