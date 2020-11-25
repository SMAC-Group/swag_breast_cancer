####################################### SWAG breast cancer ######################################################

# The aim of this script is to replicate, and make reproducible, the analysis on the breast cancer dataset of Haakensen et al. 2016

### PRELIMINARY STEPS ###

# Please refer to the readme script for the detailed passages to obtain the below object:


load("***.Rdata")


# First 31 observations in y_sub are either benign (23) or DCIS (8) and will not be considered

table(y_sub[1:31])

y <- c(rep(1,55),rep(0,70)) 

# Subset, take off first 31 lines, the X matrix wrt invasive and normal breast

complete_X <- X

X <- X[32:156,]

### NO DUPLICATES ###

# We need to use only the mean for the same name of gene to avoid duplicates.

gene <- names(table(nomen))

X_new <- matrix(rep(0,125*961),nrow = 125,ncol = 961)

for (i in 1:length(gene)) {
  
  equal <- which(nomen == gene[i])
  
  X_new[,i] <- rowMeans(X[,equal])
  
}

colnames(X_new) <- gene

# Now we take off what are not genes (all genes have 16 repetitions)

no_gene <- names(which(table(nomen) != "16"))

number <- unlist(sapply(no_gene, function(x) which(gene == x)))

colnames(X_new)[number] #quality check

# Final matrices no duplicates

X_new <- X_new[,-number]

### CV SPLIT ###

# Split: 100 for training and 25 for testing so 80% train and 20% test. 

# In particular, to preserve proportion, we need:

# train set: 44 breast cancer and 56 normal. test set: 11 breast cancer and 14 normal.

# selected seed() that guarantees preserved proportions. seed = 180

set.seed(180)

chos <- sample(index,25)  

y_test <- y[chos]

y_train <- y[-chos]

x_test <- X_new[chos,]

x_train <- X_new[-chos,]

### Data Analysis ###

# Requiring packages

dep_panning2 = c("devtools","doParallel","rngtools","doRNG","nnet","MASS","R.utils","Rcpp","RcppEigen","RcppNumerical", "panning2") 

lapply(dep_panning2, require, character.only = TRUE)


### 1st Algorithm ###

# parameters
q0 <- .05 # quantile for the whole procedure
dmax <- 8 # max model size 
mod_max <- 4e4 # model explored at each step
#nc <- 16L # number of clusters

nc <- 4

# data storage
CVs <- vector("list",dmax)
IDs <- vector("list",dmax)

VarMat <- vector("list",dmax)
set.seed(163L)
graine <- sample.int(1e6,dmax)

# Initial step: (dimension 1)
# EXAMPLE FOR D=1
cv_errors <- vector("numeric",ncol(x_train))
cl <- makeCluster(nc)
registerDoParallel(cl)

cv_errors <- foreach(i = seq_along(cv_errors), .combine = c, .packages=c("panning2")) %dopar% {
  
  seed <- graine[1] + i
  
  X <- as.matrix(cbind(rep(1,dim(x_train)[1]),x_train[,i]))
  
  y <- y_train
  
  cv <- cross_validation_logistic_count(X,y,seed)

}
stopCluster(cl)


CVs[[1]] <- cv_errors
VarMat[[1]] <- seq_along(cv_errors)

q1 <- 0.05 #quantile for first screening (see Algo1). I can be different from q0, as it selects the size of the bowl.

IDs[[1]] <- which(cv_errors <= quantile(cv_errors,q1))

id_screening <- IDs[[1]]

### General Procedure ###

# Dimension d >= 2
for(d in 2L:dmax){
  # cv0 <- cv1
  idRow <- IDs[[d-1]]
  if(d==2){
    idVar <- VarMat[[d-1]][idRow]
    nrv <- length(idVar)
  }else{
    idVar <- VarMat[[d-1]][idRow,]
    nrv <- nrow(idVar)
  }
  # build all possible 
  A <- matrix(nr=nrv*length(id_screening),nc=d)
  A[,1:(d-1)] <- kronecker(cbind(rep(1,length(id_screening))),idVar)
  A[,d] <- rep(id_screening,each=nrv)
  B <- unique(t(apply(A,1,sort)))
  id_ndup <- which(apply(B,1,anyDuplicated) == 0)
  var_mat <- B[id_ndup,]
  rm(list=c("A","B"))
  
  if(nrow(var_mat)>mod_max){
    set.seed(graine[d]+1)
    VarMat[[d]] <- var_mat[sample.int(nrow(var_mat),mod_max),]
  }else{
    VarMat[[d]] <- var_mat
  }
  
  var_mat <- VarMat[[d]]
  
  cv_errors <- rep(NA,nrow(var_mat))
  cl <- makeCluster(nc)
  registerDoParallel(cl)
  
  cv_errors <- foreach(i = seq_along(cv_errors), .combine = c, .packages=c("panning2")) %dopar% {
    
    rc <- var_mat[i,]
    cv <- NA
    seed <- graine[1] + i
    
    X <- as.matrix(cbind(rep(1,dim(x_train)[1]),x_train[,rc]))
    
    y <- y_train
    
    cv <- cross_validation_logistic_count(X,y,seed)
   
    return(cv)
  }
  stopCluster(cl)
  
  attr(cv_errors,"rng") <- NULL
  
  CVs[[d]] <- cv_errors
  cv1 <- quantile(cv_errors,probs=q0,na.rm=T)
  IDs[[d]] <- which(cv_errors<=cv1)
  
  save(CVs,file="CVs_ahus.rda")
  save(IDs,file="IDs_ahus.rda")
  save(VarMat,file="VarMat_ahus.rda")
  print(d)
}

# Graph of CV errors

m_vector <- sapply(CVs, function(x) summary(x)[4])

m_vector <- sapply(CVs, function(x) summary(x)[3]) #median

l_vector <- sapply(CVs, function(x) summary(x)[1])

u_vector <- sapply(CVs, function(x) summary(x)[6])

plot(1:length(CVs),m_vector,main = "Average CV Errors Hamming",ylab = "Average cv-error",xlab = "Model Dimension")



# Range graph of SWAG solutions

u <- 5  #upper model size (in our case at max 5 miRNAs per model) 

q0 <- 0.01

cv_target <- quantile(CVs[[u]],probs=q0,na.rm=T)

l <- ifelse(l_vector <= cv_target, 1,0) #lower model size

# Vector of dimention of the model selected 
dim_model = l:u

### Post processing ###

# Find the index of model selected and number of models in the chosen configuration

source("mod_ind_numb.R")

res <- mod_ind_numb(dim_model,CVs,cv_target)

## function to get the selected covariate name, frequency and position in the X matrix

source("list_sel_cov.R") 

# Outputs are: a list in which we have all the models of a given size and a table with all

# the names of variables present in the pool and their numerosity in the chosen models

res_2 <- list_sel_cov(VarMat,dim_model,res$sel_mod)


sel_list <- res_2$model_select

source("mat_mod.R")

mod_mat <- mat_mod(sel_list,dim_model) #matrix that contains all the 112 selected models by SWAG

source("logistic_eval.R")

output <- logistic_eval(dim_model,fil,x_train,x_test,y_train,y_test)

miRNAs_imp <- sort(table(fil)[-1]/dim(fil)[1],decreasing = T)

miRNA_prop <- as.numeric(noquote(names(miRNAs_imp)))

selection <- colnames(x_train)[miRNA_prop]

names(miRNAs_imp) <- selection

### Lasso Results ###

require(caret)

require(glmnet)

set.seed(180) # for replicability

cv.zero_mod_class <- cv.glmnet(x_train,y_train, alpha=1,family = "binomial",type.measure = "class") #a way in which you can select lambda (10k CV)

plot(cv.zero_mod_class) #it does plot the confidence intervals of the CV-errors depending on the type of measure. 

cv.zero_mod_class$cvm #mean of cv-errors (red curve)

coef_L1_min <- predict(cv.zero_mod_class,type="coefficient",s = "lambda.min")

coef_L1_1SE <- predict(cv.zero_mod_class,type="coefficient")  #default is lambda.1SE

sel <- which(coef_L1_1SE != 0) - 1 # 12 selection for 1SE excluding intercept

coef_L1_1SE[which(coef_L1_1SE != 0)]

name_las_1SE <- colnames(x_train)[sel] 

sel_min <- which(coef_L1_min != 0) - 1 # 12 selection for min excluding intercept

name_las_min <- colnames(x_train)[sel_min] #

which(name_las_1SE %in% names(var_imp)) #only name_las_min[12] not present in swag set

### Prediction test set ###

y_las_min <- as.numeric(predict(cv.zero_mod_class,type="response",s = "lambda.min",newx = x_test))

ref <- as.factor(y_test)

lasso_mat <- caret::confusionMatrix(data=y_las_min, reference = ref, positive = "1")

# Table 2 paper

lasso_resp  <- ifelse(y_las_min >= 0.5,1,0)  

lasso_resp <- factor(lasso_resp,levels = c("0","1"))

y_test <- factor(y_test,levels = c("0","1"))

conf_mat_part <- caret::confusionMatrix(data = lasso_resp,reference = y_test,positive = "1") 

# Sens. Spec. for SWAG models

test_pred <- output$fit_prob[,which(swag_auc == "1")] #54,89 best AUC

# 16/11 Finding the range for spec-sens at 0.5 cut-off

sens <- sapply(swag_roc_all_models, function(x) x[81,1]) 

spec <- sapply(swag_roc_all_models, function(x) 1 - x[81,2])

fcm <- output$conf_mat

fcm[[1]]$byClass

ppv <-  sapply(fcm, function(x) x$byClass[3])

npv <- sapply(fcm, function(x) x$byClass[4])

acc <- sapply(fcm, function(x) x$overall[1])

swag_mat <- cbind(acc,sens, spec, npv, ppv)

colnames(swag_mat) <- c("Accuracy", "Sensitivity","Specificity","Negative Predicted Value", "Positive Predicted Value")

rownames(swag_mat) <- c()


