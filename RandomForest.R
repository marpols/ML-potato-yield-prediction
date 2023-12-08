#****************************************************************************
#Random Forest on field data provided by Matt Ramsay, Climate data from ECCC
#and soil data from field and CANSIS
#Created by Mariaelisa Polsinelli for AAFC, 2023
#****************************************************************************

if (!require("randomForest")){
  install.packages("randomForest")
}
if (!require("caret")){
  install.packages("caret")
}
if (!require("dplyr")){
  install.packages("dplyr")
}
if (!require("Metrics")){
  install.packages("Metrics")
}
if (!require("ggplot2")){
  install.packages("ggplot2")
}
packages <- c("dplyr","readxl","randomForest","caret","Metrics","ggplot2",
              "tidyr","reshape2","RColorBrewer")
lapply(packages, library, character.only=TRUE)

source("rfFunctions.R")
source("rfPlots.R")

dataname <- "MR"
run <- "noNormcv5"

kfolds <- 5 #number of folds
nsize <- c(1,5,10) #min nodesize -> larger number=smaller tree. Regression default=5

seed <- 1234

hybrid <- F
normalization <- F
normtype <- 1 #1 for log scale, 2 for min/max

if (hybrid){
  load()
  dir <- paste(wd,"/",dataname,"_Hybrid",run,sep="")
  dir.create(dir)
  setwd(dir)
} else{
  load("MLalone_original.RData")
  dir <- paste(wd,"/",dataname,"_RF",run,sep="")
  dir.create(dir)
  setwd(dir)
}

#RANDOM FOREST------------------------------------------------------------------

i <- 1
begin <- Sys.time()

for (f in fields){
  dir.create(f)
  sink(file=paste(f,"/",f,"_",run,"_log.txt",sep=""), split = TRUE) #write output to log
  set.seed(seed)
  print(paste("running rf for",f,"at",Sys.time(),sep=" "))
  
  #Separate each individual field year out as test data
  ss_data <-  ss_data[sample(1:nrow(ss_data)), ] 
  train_x <- ss_data[ss_data$CFIDYr != f, 2:c]
  all_train <- ss_data[ss_data$CFIDYr != f, 1:c] %>% select(1,4:c)
  train_y <- ss_data[ss_data$CFIDYr != f, 1]
  test_x <- ss_data[ss_data$CFIDYr == f, 2:c]
  test_y <- ss_data[ss_data$CFIDYr == f,1]
  
  if (normalization){
    if (normtype==1){
      train_x <- log(as.data.frame(train_x))
      all_train <- log(as.data.frame(all_train))
      test_x <- log(as.data.frame(train_x))
    }else{
      process_trn <- preProcess(as.data.frame(train_x), method=c("range"))
      process_alltrn <- preProcess(as.data.frame(all_train), method=c("range"))
      process_tst <- preProcess(as.data.frame(test_x), method=c("range"))
      train_x <- predict(process, as.data.frame(train_x))
      all_train <- predict(process, as.data.frame(all_train))
      test_x <- predict(process, as.data.frame(test_x))
    }
  }
  
  train_yrs <- train_x[,2]
  test_yrs <- test_x[,2]
  train_flds <- train_x[,1]
  test_flds <- test_x[,1]
  
  train_x <- train_x[,3:length(train_x)]
  test_x <- test_x[,3:length(test_x)]
  
  #feature selection and hyperparameter tuning
  features <- featselec(train_x,train_y,test_x,test_y,fold=kfolds)
  fnum <- min_mse(features) #get optimal number of features
  write.csv(fnum$features,paste(f,"/features_",run,".csv",sep=""))
  
  #update train and test sets to have optimum number of features
  train_x <- train_x[which(colnames(train_x) %in% fnum$features)]
  test_x <- test_x[which(colnames(test_x) %in% fnum$features)]
  all_train <- all_train[which(colnames(all_train) %in% c("yield_tha",fnum$features))]
  
  #run rf with grid search, cv, to find best model
  s <- c(1:fnum$metrics$num_features)
  if(is.null(nsize)){ #grid search for mtry and ntree
    rfs <- hype_tune(all_train,test_x,test_y,folds=kfolds,mtry_seq=s)
    best <- best_mod(rfs)
    prediction <- best$predictions
    train_pred <- best$train_pred
  } else{ #gridserch for mtry, ntree and nodesize
    min_nodesrch <- list()
    for (minn in nsize){
      print(paste("------ Nodesize: ",minn,sep=""))
      rfs <- hype_tune(all_train,test_x,test_y,folds=kfolds,
                       mtry_seq=s, nodesize = minn)
      best <- best_mod(rfs)
      min_nodesrch[[toString(minn)]] <- best
    }
    best <- best_mod(min_nodesrch,nodesize=T)
    prediction <- best$model$predictions
    train_pred <- best$model$train_pred
  }
  
  #save all predicted vs. observed field-year yield data points as .csv
  field_pred <- data.frame(labels[i],test_y,prediction)
  colnames(field_pred) <- c("Label","Actual","Predicted")
  if (i == 1){
    write.table(field_pred, file=paste("allFields_",run,".csv",sep=""), 
                row.names = F,col.names = T, append = F, sep=",")
  } else {
    write.table(field_pred, file=paste("allFields_",run,".csv",sep=""), 
                row.names = F,col.names = F, append = T, sep=",")
  }
  
  write.csv(field_pred, file=paste(dir,"/",f,"/preds_",run,".csv",sep=""))
  
  #save metrics
  train_metrics <- calc_metrics(train_y,train_pred)
  test_metrics <- calc_metrics(test_y, prediction)
  
  metrics <- data.frame(Training=train_metrics,
                        Testing=test_metrics)
  metrics <- rbind(metrics,c(c(".", best$ntree)))
  metrics <- rbind(metrics,c(c(".", best$mtry)))
  metrics <- rbind(metrics,c(c(".", length(prediction))))
  rownames(metrics) <- c("MSE","RMSE","RRMSE","R2","BIAS","PBIAS",
                         "ntree","mtry","size.test.data")
  write.csv(metrics,paste(dir,"/",f,"/metrics.csv",sep=""))
  
  #get average predicted and observed yield for field-year
  avg_pred <- mean(prediction)
  avg_obs <- mean(test_y)
  
  overall <- data.frame(CFIDYr=labels[i],Year=test_yrs[1],
                        Actual=avg_obs,Predicted=avg_pred)
  if (i == 1){
    write.table(overall, file=paste("ActualVsPredicted_avg",run,".csv",sep=""), 
                row.names = F,col.names = T, append = F, sep=",")
  } else {
    write.table(overall, file=paste("ActualVsPredicted_avg",run,".csv",sep=""), 
                row.names = F,col.names = F, append = T, sep=",")
  }
  i <- i + 1
  print(paste("Ended at",Sys.time(),sep=" "))
  sink()
}
print(Sys.time() - begin)

generate_plot(paste("ActualVsPredicted_avg",run,".csv",sep=""),run,dataname)