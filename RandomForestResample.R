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
if (!require("UBL")){
  install.packages("UBL")
}
packages <- c("dplyr","readxl","randomForest","caret","Metrics","ggplot2",
              "tidyr","reshape2","RColorBrewer","beepr","UBL")
lapply(packages, library, character.only=TRUE)

source("rfFunctions.R")
source("rfPlots.R")

dataname <- "MR"
run <- "SMOTER_B-cv2"
ptitle <- "Random Forest Balanced OverSampling" #title for act v. pred plot

kfolds <- 2 #number of folds
nsize <- c(1,5,10) #min nodesize -> larger number=smaller tree. 
                   #Regression default=5. Set as NULL to bypass nodesize gridsearch

seed <- 1234

hybrid <- F
normalization <- F


if(normalization){
  run <- paste("wNorm_",run,sep="")
}

if (hybrid){
  load("hybrid_lidar.RData")
  ss_data <- select(ss_data,1:107,144:145)
  c <- length(ss_data)
  dir <- paste(dataname,"_Hybrid",run,sep="")
} else{
  load("MLalone_slope-lidar2.RData")
  dir <- paste(dataname,"_RF",run,sep="")
}

dir.create(dir)
setwd(dir)

#ss_data <- subset(ss_data, select=-c(5,6)) #remove both GDDcum and GDDcum_gs
#ss_data <- subset(ss_data, select=-c(6)) #remove GDDcum_gs
ss_data <- subset(ss_data, select=-c(5:12)) #remove cumulative climate features
#ss_data <- subset(ss_data, select=-c(5,9:28,37:39,83)) #keep only GS cumulative climate features
ss_data <- ss_data %>% select(yield_tha,CFIDYr, plant_ye,everything())
#ss_data <- subset(ss_data,select=-c(6:13,42:82,85:110))
#ss_data <- subset(ss_data,select=-c(6:19,42:82,103:108)) #remove GDD vars, keep only min and max temp vars up to october

c <- length(ss_data)

#RANDOM FOREST------------------------------------------------------------------
# run '#######' for specific individual field-years (not using for loop) 
i <- 1

#i <- match("OffeyCornerRight-2019",fields) #######
f <- fields[i] #######
#fname <- paste(f,"_2",sep="")

begin <- Sys.time()

#for (f in fields){
  dir.create(f)
  sink(file=paste(f,"/",f,"_",run,"_log.txt",sep=""), split = TRUE) #write output to log
  
  print(paste("running rf for",f,"at",Sys.time(),sep=" "))
  
  #Separate individual field year out as test data
  set.seed(seed)
  all_train <- ss_data[ss_data$CFIDYr != f,] %>% select(1,4:c) #training dataset including Y
  
  
  #Resampling of training data
  type <- "SmoteRegress" #calls = RandUnderRegress, RandOverRegress, SmoteRegress
  x <- all_train
  
  rpoints <- boxplot.stats(x[,1])
  rpoints$stats
  #originial distribution plots
  ggplot(x) + geom_boxplot(aes(x=yield_tha))
  ggsave(paste(f,"/org_boxplot.png",sep=""))
  ggplot(x) + geom_histogram(aes(x=yield_tha),bins=26)
  ggsave(paste(f,"/org_hist.png",sep=""))
  
  # #3 points
  # rel <- t(matrix(data=c(c(rpoints$stats[1],1.0,0.0),
  #                        c(rpoints$stats[3],0.0,0.0),
  #                        c(rpoints$stats[5],1.0,0.0)),
  #                 nrow=3,ncol=3))
  
  # #4 points
  # rel <- t(matrix(data=c(c(rpoints$stats[1],1.0,0.0),
  #                        c(39,0.4,0.0),
  #                        c(rpoints$stats[3],0.0,0.0),
  #                        c(rpoints$stats[5],0.0,0.0)),
  #                 nrow=3,ncol=4))
  
  #5 points
  rel <- t(matrix(data=c(c(rpoints$stats[1],1.0,0.0),
                         c(40,0.4,0.0),
                         c(rpoints$stats[3],0.0,0.0),
                         c(50,1.0,0.0),
                         c(55,0.0,0.0)),
                  nrow=3,ncol=5))

  C.perc = list(3.5,0.7,1.0,0.7) #UNder: <= 1, Over: =>1, Smote: (Over,Under) or (Over,Under,OVer) 
  k <- NA   #Nearest neighbours for SMOTER
  
  set.seed(seed)
  new_data <- if(!is.na(k)) 
    do.call(type,list(yield_tha~.,x,C.perc=C.perc,rel,k)) else do.call(type,list(yield_tha~.,x,C.perc=C.perc,rel))
  
  #new distributions
  ggplot(new_data) + geom_boxplot(aes(x=yield_tha))
  ggsave(paste(f,"/",type,"_boxplot.png",sep=""))
  ggplot(new_data) + geom_density(aes(x=yield_tha))
  ggsave(paste(f,"/",type,"_density.png",sep=""))
  ggplot(new_data) + geom_histogram(aes(x=yield_tha),bins=26)
  ggsave(paste(f,"/",type,"_hist.png",sep=""))
  
  new_points <- boxplot.stats(new_data[,1])
  new_points$stats
  
  #save resampled data and parameters
  write.csv(list(org.quantiles = rpoints$stats,
                 new.quantiles=new_points$stats),
            paste(f,"/resamplingQuant_",run,".csv",sep=""))
  parameters <- rbind(rel,c(C.perc,'.','.'),c(k,'.','.'))
  
  # #3 points
  # rownames(parameters) <- c('rel.min','rel.med','rel.max','C.perc','k (SMOTER)')
  # #4 points
  # rownames(parameters) <- c('rel.min','rel.2q','rel.med','rel.max','C.perc','k (SMOTER)')
  # #5 points
  rownames(parameters) <- c('rel.min','rel.2q','rel.med','rel.3q','rel.max','C.perc','k (SMOTER)')
  
  write.csv(parameters,
            paste(f,"/parameters_",run,".csv",sep=""))
  write.csv(new_data,
            paste(f,"/newtrain_",run,".csv",sep=""))

  all_train <- new_data
  #separate train and test x and y
  train_x <- all_train[,2:length(all_train)]
  train_y <- all_train[, 1]
  test_x <- ss_data[ss_data$CFIDYr == f, 2:c]
  test_y <- ss_data[ss_data$CFIDYr == f,1]
  
  #set aside year and field names
  test_yrs <- test_x[,2]
  test_flds <- test_x[,1]
  
  #final test sets
  test_x <- test_x[,3:length(test_x)]
  
  if (normalization){
    train_x <- norm_train(train_x)
    all_train <- norm_alltrain(train_x, train_y)
    test_x <- norm_test(test_x)
  }
  
  #f <- paste(f,"_2",sep="") #######
  
  #feature selection and hyperparameter tuning
  features <- featselec(train_x,train_y,test_x,test_y,fold=kfolds,af=10)
  beep(sound=4)
 
  #fnum <- min_mse(features) #get optimal number of features - need to update function
  # set fnum <- features[[num features]] to choose manually #######
  
  write.csv(list(features=fnum$features,per_increaseMSE=fnum$MSE_importance,
                 NodePurity=fnum$NodeImp_importance),
                 paste(f,"/features_",run,".csv",sep=""))
  print(paste("Optimum number of features:", fnum, sep=" "))
  
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
  } else{ #gridsearch for mtry, ntree and nodesize
    min_nodesrch <- list()
    for (minn in nsize){
      print(paste("------ Nodesize: ",minn,sep=""))
      rfs <- hype_tune(all_train,test_x,test_y,folds=kfolds,
                       mtry_seq=s, nodesize=minn)
      best <- best_mod(rfs)
      min_nodesrch[[toString(minn)]] <- best
    }
    best <- best_mod(min_nodesrch,nodesize=T)
    prediction <- best$model$predictions
    train_pred <- best$model$train_pred
  }
  beep(sound=4)
  
  #save all predicted vs. observed field-year yield data points as .csv
  field_pred <- data.frame(labels[i],test_y,prediction)
  colnames(field_pred) <- c("Label","Actual","Predicted")
  
  colnames <- if(i == 1) T else F
  append <- if(i == 1) F else T 
 
  write.table(field_pred, file=paste("allFields_",run,".csv",sep=""), 
              row.names = F,col.names = colnames, append = append, sep=",")

  
  write.csv(field_pred, file=paste(getwd(),"/",f,"/preds_",run,".csv",sep=""))
  
  #save metrics
  train_metrics <- calc_metrics(train_y,train_pred)
  test_metrics <- calc_metrics(test_y, prediction)
  
  metrics <- data.frame(Training=train_metrics,
                        Testing=test_metrics)
  if(is.null(nsize)){
    metrics <- rbind(metrics,c(c(".", best$ntree)))
    metrics <- rbind(metrics,c(c(".", best$mtry)))
    metrics <- rbind(metrics,c(c(".", length(fnum$features))))
    metrics <- rbind(metrics,c(c(".", length(prediction))))
    rownames(metrics) <- c("MSE","RMSE","RRMSE","R2","BIAS","PBIAS",
                           "ntree","mtry","num.features","size.test.data")
  }else{
    metrics <- rbind(metrics,c(c(".", best$nodesize))) 
    metrics <- rbind(metrics,c(c(".", best$model$ntree)))
    metrics <- rbind(metrics,c(c(".", best$model$mtry)))
    metrics <- rbind(metrics,c(c(".", length(fnum$features))))
    metrics <- rbind(metrics,c(c(".", length(prediction))))
    rownames(metrics) <- c("MSE","RMSE","RRMSE","R2","BIAS","PBIAS","nodesize",
                           "ntree","mtry","num.features","size.test.data")
  }
  write.csv(metrics,paste(getwd(),"/",f,"/metrics.csv",sep=""))
  
  #get average predicted and observed yield for field-year
  avg_pred <- mean(prediction)
  avg_obs <- mean(test_y)
  
  overall <- data.frame(CFIDYr=labels[i],Year=test_yrs[1],
                        Actual=avg_obs,Predicted=avg_pred)

  write.table(overall, file=paste("ActualVsPredicted_avg",run,".csv",sep=""), 
              row.names = F,col.names = colnames, append = append, sep=",")
  
  i <- i + 1
  print(paste("Ended at",Sys.time(),sep=" "))
  print(Sys.time() - begin)
  sink()
#}
  
print(Sys.time() - begin)

generate_plot(paste("ActualVsPredicted_avg",run,".csv",sep=""),run,dataname,title=ptitle)

beep(sound=8)
beep(sound=4)
