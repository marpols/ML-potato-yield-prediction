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
packages <- c("dplyr","readxl","randomForest","caret","Metrics","ggplot2","tidyr","reshape2","RColorBrewer")
lapply(packages, library, character.only=TRUE)

wd <- getwd()
mlalone <- paste(wd,"/MattRamsay_RU_20152021/mrData.csv",sep="")
#hybrid_data <- paste(wd,"",sep="")
dataname <- "MR"
run <- "noNormcv5"

seed <- 1234

hybrid <- F
normalization <- F
normtype <- 1 #1 for log scale, 2 for min/max

if (hybrid){
  data <- data.frame(read.csv(hybrid_data))
  dir <- paste(wd,"/",dataname,"_Hybrid",run,sep="")
  dir.create(dir)
  setwd(dir)
} else{
  data <- data.frame(read.csv(mlalone))
  dir <- paste(wd,"/",dataname,"_RF",run,sep="")
  dir.create(dir)
  setwd(dir)
}

#select only datapoints that have %sand info (have soil samples)
ss_data <- subset(data, data$sh_PER_SAN != 0) %>% select(6:33,35,41:47,50:55,58,59,62:69,72:79,82:89,92:108)
# yld_dist_snd <- ggplot(ss_data, aes(x=yield_tha)) + geom_histogram(color="black", fill="white") + labs(x="Yield (t/ha)", title = "Distribution of Russet Yield Datapoints with Field Soil Information")
# yld_dist_snd
ss_data$PH1[which(ss_data$PH1 %in% c(0))] <- 6.0 #replace missing pH values with pH of 6 as per CANSIS data

#move field ID, year, and yield columns to front
ss_data <- ss_data %>% dplyr::select("plant_ye", everything())
ss_data <- ss_data %>% dplyr::select("CFIDYr", everything())
ss_data <- ss_data %>% dplyr::select("yield_tha", everything()) 

#functions----------------------------------------------------------------------
calc_metrics <- function(actual,predicted){
  mse <- mse(actual, predicted)
  rmse <- rmse(actual, predicted)
  rrmse <- rmse/mean(actual)
  r2 <- cor(actual,predicted)^2*100
  bias <- bias(actual, predicted)
  pbias <- bias/mean(actual)*100
  return (c(mse,rmse,rrmse,r2,bias,pbias))
}

min_mse <- function(fs_results){
  
  min_train <- fs_results[[1]]$metrics$train_mse
  train_n <- 1
  min_test <- fs_results[[1]]$metrics$test_mse
  test_n <- 1
  
  i <- 2
  while (i <= length(fs_results)){
    if (fs_results[[i]]$metrics$train_mse < min_train){
      min_train <- fs_results[[i]]$metrics$train_mse
      train_n <- i
    }
    if (fs_results[[i]]$metrics$test_mse < min_test){
      min_test <- fs_results[[i]]$metrics$test_mse
      test_n <- i
    }
    i <- i + 1
  }
  min_errors <- list(train_mse=min_train,train_vars=train_n,
               test_mse=min_test,test_vars=test_n)
  a <- features[[min_errors$test_vars]]$metrics$train_mse
  b <- features[[min_errors$train_vars]]$metrics$test_mse

  diff_test <- min_errors$test_mse - a
  diff_train <- b - min_errors$train_mse
  
  if(diff_test < diff_train & diff_test > 0){
    return (fs_results[[test_n]])}else{return (fs_results[[train_n]])}
}

best_mod <- function(models){
  btest <- 1000 #arbitrary initial value
  bntree <- 100 
  bmtry <- 10
  bpreds <- NULL
  bmod <- NULL
  for (mod in models){
    if  (mod$test_metrics[2] < btest){
      btest <- mod$test_metrics[2] 
      bntree <- mod$rf$dots$ntree
      bmtry <- mod$rf$bestTune$mtry
      bpreds <- mod$prediction
      bmod <- mod
      }
  }
  return (list(ntree=bntree,mtry=bmtry,train_pred=bmod$rf$finalModel$predicted,
               predictions=bpreds,
               test_metrics=list(mse=bmod$test_metrics[1],
                                 rmse=bmod$test_metrics[2],
                                 rrmse=bmod$test_metrics[3],
                                 r2=bmod$test_metrics[4],
                                 bias=bmod$test_metrics[5],
                                 pbias=bmod$test_metrics[6])))
}

#Feature Selection--------------------------------------------------------------
featselec <- function(train_x,train_y,test_x,test_y, nt=500, fold=10){
  
  set.seed(seed)
  
  varlist <- list()
  
  #CV for feature selection
  print("getting optimum number of features")
  CV <- rfcv(
    train_x,
    train_y,
    cv.fold = fold,
    scale = "log",
    step = 0.5,
    mtry = function(p)
      max(1, floor(sqrt(p))))
  
  nvars <- as.data.frame(CV$error.cv)
  nvars$vars <- as.numeric(rownames(nvars))
  v <- nvars$vars[nvars$`CV$error.cv` == min(nvars$`CV$error.cv`)]
  
  mtry <- as.data.frame(tuneRF(train_x, train_y, stepFactor = 1.5, trace=F, plot=F))
  n <- mtry$mtry[mtry$OOBError == min(mtry$OOBError)]
  
  print("running random forest")
  rf <- randomForest(
    train_x,
    y = train_y,
    proximity = TRUE,
    ntree = nt,
    mtry = n,
    plot = TRUE)
  print(rf)
  
  print("getting variable importance")
  #get variable importance
  var_imp <- rf$importance
  vars <- row.names(rf$importance)
  importance <- data.frame(vars, var_imp) %>% arrange(desc(var_imp))
  
  print("starting add-one-in feature selection approach")
  j <- 1
  while (j <= v + 10) { #try additional 10 features
    set.seed(seed)
    trnx <- train_x[importance$vars[1:j]]
    tstx <- test_x[importance$vars[1:j]]
    mtry <- as.data.frame(tuneRF(trnx, train_y, stepFactor = 1.5))
    n <- mtry$mtry[mtry$OOBError == min(mtry$OOBError)]
    rf <-
      randomForest(
        trnx,
        y = train_y,
        proximity = TRUE,
        ntree = nt,
        mtry = n,
        plot = TRUE
      )
    trnm <- calc_metrics(train_y, rf[3]$predicted)
    pred <- predict(rf, tstx, type = "response", predict.all = TRUE)
    tstm <- calc_metrics(test_y, pred[[1]])
    
    #number of features, ntree, mtry, train mse, test mse
    val <- list("num_features"=j, "ntrees"=nt, "mtry"=n,"train_mse"=trnm[1],"test_mse"=tstm[1])
    info <- list(metrics=val,features=importance$vars[1:j])
    varlist[[j]] <- info

    print(paste("number of features:",j,"ntree:",nt,"mtry:",n,"train MSE:",
                trnm[1],"test MSE:",tstm[1],sep=" "))
    print(importance$vars[1:j])
    
    j <- j + 1
  }
  return (varlist)
}

#Hyperparameter tuning (mtry and ntree)-----------------------------------------
hype_tune <- function(train, test, y, folds=10, repeats=3, mtry_seq=(1:60), 
                     ntree_seq=seq(100,1000,by=100),tune_search='grid'){
#train - all training data (including train y)
#test - test X values 
#y - test y values
  starttime <- Sys.time()
  modellist <- list()
  
  control <- trainControl(method='repeatedcv', 
                          number=folds, 
                          repeats=3, 
                          search=tune_search)
  
  tunegrid <- expand.grid(.mtry = mtry_seq)
  
  print(paste("Training with ",folds,"-Fold Cross Validation",sep=""))
  for (nt in ntree_seq){
    set.seed(seed)
    key <- toString(nt)
    print(paste("Ntree: ",key,sep=""))
    
    rf_gridsearch <- train(yield_tha ~ ., 
                           data = train,
                           method = 'rf',
                           metric = 'RMSE',
                           trControl = control,
                           tuneGrid = tunegrid,
                           ntree=nt)
    
    rf_pred <- predict(rf_gridsearch,test)
    pred_res <- calc_metrics(y,as.numeric(rf_pred))
    
    n <- as.numeric(rf_gridsearch$bestTune$mtry)
    print(paste("best mtry:",n, 
                "train_RMSE:",rf_gridsearch$results$RMSE[rf_gridsearch$results$mtry==n], 
                "test_RMSE:", pred_res[2],sep=" "))
    print(Sys.time() - starttime)
    info <- list(rf = rf_gridsearch, prediction=rf_pred, test_metrics=pred_res)
    modellist[[key]] <- info
    
  }
  return (modellist)
}

#RANDOM FOREST------------------------------------------------------------------

c <- length(ss_data)

fields <- unique(ss_data$CFIDYr) %>% sort() #get list of full field names
labels <- c("L-2015","OS-2017","I17-2016","L-2017","CN-2017","P-2017","BR-2017",
            "JM-2017","M1-2017","I13-2017","OCR-2019","IC-2019","I17-2019",
            "RB-2019","P-2019","M4-2017","OS-2019","L-2021","JM-2020",
            "GRS-2019") %>% sort()

#Seperate each individual field year out as test data
i <- 1
begin <- Sys.time()
for (f in fields){
  set.seed(seed)
  dir.create(f)
  print(paste("running rf for",f,"at",Sys.time(),sep=" "))
  
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
  features <- featselec(train_x,train_y,test_x,test_y)
  fnum <- min_mse(features) #get optimal number of features
  write.csv(fnum$features,paste(f,"/features_",run,".csv",sep=""))
  
  #update train and test sets to have optimum number of features
  train_x <- train_x[which(colnames(train_x) %in% fnum$features)]
  test_x <- test_x[which(colnames(test_x) %in% fnum$features)]
  all_train <- all_train[which(colnames(all_train) %in% c("yield_tha",fnum$features))]
  
  #run rf with grid search, cv, to find best model
  s <- fnum$metrics$num_features
  rfs <- hype_tune(all_train,test_x,test_y,folds=5,mtry_seq=(1:s))
  
  best <- best_mod(rfs)
  prediction <- best$predictions
  
  #save all predicted vs. observed field-year yield data points as .csv
  field_pred <- data.frame(labels[i],test_y,prediction)
  colnames(field_pred) <- c("Label","Actual","Predicted")
  if (i == 1){
    write.table(field_pred, file=paste("allFields_",run,".csv",sep=""), row.names = F,col.names = T, append = F, sep=",")
  } else {
    write.table(field_pred, file=paste("allFields_",run,".csv",sep=""), row.names = F,col.names = F, append = T, sep=",")
  }
  
  write.csv(field_pred, file=paste(dir,"/",f,"/preds_",run,".csv",sep=""))
  
  train_metrics <- calc_metrics(train_y,best$train_pred)
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
    write.table(overall, file=paste("ActualVsPredicted_avg",run,".csv",sep=""), row.names = F,col.names = T, append = F, sep=",")
  } else {
    write.table(overall, file=paste("ActualVsPredicted_avg",run,".csv",sep=""), row.names = F,col.names = F, append = T, sep=",")
  }
  i <- i + 1
  print(paste("Ended at",Sys.time(),sep=" "))
}
print(Sys.time() - begin)

#Overall Metrics and Plot------------------------------------------------------
avg_fields <- read.csv(paste("ActualVsPredicted_avg",run,".csv",sep=""))

metrics_oa <- calc_metrics(avg_fields$Actual, avg_fields$Predicted)

metrics_all <- data.frame(metrics_oa)
colnames(metrics_all) <- c("Overall Field-Years")
rownames(metrics_all) <- c("MSE","RMSE","RRMSE","R2","BIAS","PBIAS")
write.csv(paste(dataname,"_metrics",run,".csv",sep=""))

#generate plot
plot <- ggplot(avg_fields) + geom_point(aes(x=Actual, y = Predicted, 
                                           color=factor(Year), size=1),shape=1) + geom_abline(intercept=0, slope=1) +
  labs(x = "Actual yield (t/ha.)",y="Predicted yield (t/ha.)", title="Random Forest Actual Vs. Predicted Yield", color="Year") +
  scale_x_continuous(breaks = seq(0,60, by=5)) + scale_y_continuous(breaks = seq(0,60, by=5)) + coord_cartesian(ylim = c(27,55), xlim=c(27,55)) +
  theme(text = element_text(size = 18),legend.text=element_text(size=10)) +
  guides(size="none",colour = guide_legend(override.aes = list(size=5)))

ggsave(paste(dataname,"_actVpredPlot",run,".png",sep=""), width=7, height=5)

