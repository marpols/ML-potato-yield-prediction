#****************************************************************************
#Random Forest on field data provided by Matt Ramsay, Climate data from ECCC
# and soil data from field and CANSIS
# Created by Mariaelisa Polsinelli for AAFC, 2023
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
run <- 1


hybrid <- F

if (hybrid){
  data <- data.frame(read.csv(hybrid_data))
  dir <- paste(wd,"/",dataname,"_Hybridrun",run,sep="")
} else{
  data <- data.frame(read.csv(mlalone))
  dir <- paste(wd,"/",dataname,"_RFrun",run,sep="")
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

#Feature Selection--------------------------------------------------------------
featselec <- function(train_x,train_y,test_x,test_y, nt=500){
  
  set.seed(1234)
  
  varlist <- list()
  
  #CV for feature selection
  CV <- rfcv(
    train_x,
    train_y,
    cv.fold = 5,
    scale = "log",
    step = 0.5,
    mtry = function(p)
      max(1, floor(sqrt(p))))
  
  nvars <- as.data.frame(CV$error.cv)
  nvars$vars <- as.numeric(rownames(nvars))
  v <- nvars$vars[nvars$`CV$error.cv` == min(nvars$`CV$error.cv`)]
  
  mtry <- as.data.frame(tuneRF(train_x, train_y, stepFactor = 1.5))
  n <- mtry$mtry[mtry$OOBError == min(mtry$OOBError)]
  
  rf <- randomForest(
    train_x,
    y = train_y,
    proximity = TRUE,
    ntree = nt,
    mtry = n,
    plot = TRUE)
  print(rf)
  
  #get variable importance
  var_imp <- rf$importance
  vars <- row.names(rf$importance)
  importance <- data.frame(vars, var_imp) %>% arrange(desc(var_imp))
  write.csv(importance, file = "variable_importance.csv", row.names = F)
  
  j <- 1
  while (j <= v + 5) { #use additional 5 features
    trnx <- train_x[importance$vars[1:j]]
    tstx <- test_x[importance$vars[1:j]]
    mtry <- as.data.frame(tuneRF(train_x, train_y, stepFactor = 1.5))
    n <- mtry$mtry[mtry$OOBError == min(mtry$OOBError)]
    rf <-
      randomForest(
        trnx,
        y = train_y,
        proximity = TRUE,
        ntree = 300,
        mtry = n,
        plot = TRUE
      )
    trnm <- calc_metrics(train_y, rf[3]$predicted)
    pred <- predict(rf, tstx, type = "response", predict.all = TRUE)
    tstm <- calc_metrics(test_y, pred[[1]])
    
    #number of features, ntree, mtry, train mse, test mse
    val <- data.frame(v, nt, n, trnm[1], testm[1])
    varlist[[j]] <- val
    
    print(paste("number of features:",v,"ntree:",nt,"mtry:",n,"train MSE:",
                trn[1],"test MSE:",testm[1],sep=" "))
    print(importance$vars[1:j])
    
    j <- j + 1
  }
  return (varlist)
}


#Hyperparameter tuning (mtry and ntree)-----------------------------------------
hype_tune <- function(train, test, y, folds=10, repeats=3, mtry_seq=(1:60), 
                     ntree_seq=seq(100,900,by=200)){
#train - all training data (including train y)
#test - test X values 
#y - test y values
  
  modellist <- list()
  
  control <- trainControl(method='repeatedcv', 
                          number=folds, 
                          repeats=3, 
                          search='grid')
  
  tunegrid <- expand.grid(.mtry = mtry_seq)
  
  print(paste("Training with ",folds,"-Fold Cross Validation",sep=""))
  for (nt in ntree_seq){
    set.seed(1234)
    key <- toString(nt)
    print(paste("Ntree: ",key,sep=""))
    
    rf_gridsearch <- train(yield_tha ~ ., 
                           data = train,
                           method = 'rf',
                           metric = 'RMSE',
                           trControl = control,
                           tuneGrid = tunegrid,
                           ntree=nt)
    
    modellist[[key]] <- rf_gridsearch
    
    rf_pred <- predict(rf_gridsearch,test)
    pred_res <- calc_metrics(y,as.numeric(rf_pred))
    
    n <- as.numeric(rf_gridsearch$bestTune$mtry)
    print(paste("best mtry:",n, 
                "train_RMSE:",rf_gridsearch$results$RMSE[rf_gridsearch$results$mtry==n], 
                "test_RMSE:", pred_res[2],sep=" "))
    
  }
  return (modellist)
}



#RANDOM FOREST------------------------------------------------------------------

c <- length(ss_data)

fields <- unique(ss_data$CFIDYr) %>% sort() #get list of full field names
labels <- c("L-2015","OS-2017","I17-2016","L-2017","CN-2017","P-2017","BR-2017",
            "JM-2017","M1-2017","I13-2017","OCR-2019","IC-2019","I17-2019",
            "RB-2019","P-2019","M4-2017","OS-2019","L-2021","JM-2020","GRS-2019") %>% sort()

#Seperate each individual field year out as test data
i <- 1
for (f in fields){
  train_x <- ss_data[ss_data$CFIDYr != f, 2:c]
  all_train <- ss_data[ss_data$CFIDYr != f, 1:c] %>% select(1,4:c)
  train_y <- ss_data[ss_data$CFIDYr != f, 1]
  test_x <- ss_data[ss_data$CFIDYr == f, 2:c]
  test_y <- ss_data[ss_data$CFIDYr == f,1]
  
  train_yrs <- train_x[,2]
  test_yrs <- test_x[,2]
  train_flds <- train_x[,1]
  test_flds <- test_x[,1]
  
  train_x <- train_x[,3:length(train_x)]
  test_x <- test_x[,3:length(test_x)]
  
  set.seed(1234)
  #feature selection and hyperparameter tuning
  features <- featselec(train_x,test_x,train_y,test_y)
  params <- hype_tune(all_train,test_x,test_y)
  
  # mtry <- as.data.frame(tuneRF(train_x,train_y, stepFactor=1.5))
  # n <- mtry$mtry[mtry$OOBError == min(mtry$OOBError)]
  
  #run random forest (training)
  rf <- randomForest(train_x,y=train_y, proximity=TRUE, ntree=300, mtry=n, plot=TRUE)
  
  print(rf)
  
  train_metrics <- calc_metrics(train_y,rf[3]$predicted)
  
  #Testing (predict on test set)
  prediction <- predict(rf,test_x,type="response",predict.all=TRUE)
  
  #save all predicted vs. observed field-year yield data points as .csv
  field_pred <- data.frame(labels[i],test_y,prediction[[1]])
  colnames(field_pred) <- c("Label","Actual","Predicted")
  if (i == 1){
    write.table(field_pred, file=paste(dir,"/","allFields_",run,".csv",sep=""), row.names = F,col.names = T, append = F, sep=",")
  } else {
    write.table(field_pred, file=paste(dir,"/","allFields_",run,".csv",sep=""), row.names = F,col.names = F, append = T, sep=",")
  }
  
  write.table(field_pred, file=paste(dir,"/",f,run,".csv",sep=""), row.names = F,col.names = T, append = F, sep=",")
  
  test_metrics <- calc_metrics(test_y, prediction[[1]])
  
  metrics <- data.frame(train_metrics,test_metrics)
  colnames(metrics) <- c("Training","Testing")
  rownames(metrics) <- c("MSE","RMSE","RRMSE","R2","BIAS","PBIAS")
  write.csv(paste(dir,"/",f,"metrics.csv",sep=""))
  
  #get average predicted and observed yield for field-year
  avg_pred <- mean(prediction[[1]])
  avg_obs <- mean(test_y)
  
  overall <- data.frame(labels[i],test_yrs[1],avg_obs,avg_pred)
  colnames(overall) <- c("CFIDYr","Year","Actual","Predicted")
  if (i == 1){
    write.table(overall, file=paste(dir,"/","ActualVsPredicted",run,".csv",sep=""), row.names = F,col.names = T, append = F, sep=",")
  } else {
    write.table(field_pred, file=paste(dir,"/","ActualVsPredicted_avg",run,".csv",sep=""), row.names = F,col.names = F, append = T, sep=",")
  }
  i <- i + 1
}

#Overall Metrics and Plot------------------------------------------------------
avg_fields <- read.csv(paste(dir,"/","ActualVsPredicted",run,".csv",sep=""))

metrics_oa <- calc_metrics(avg_fields$Actual, avg_fields$Predicted)

metrics_all <- data.frame(metrics_oa)
colnames(metrics_all) <- c("Overall Field-Years")
rownames(metrics_all) <- c("MSE","RMSE","RRMSE","R2","BIAS","PBIAS")
write.csv(paste(dir,"/",dataname,"_metrics",run,".csv",sep=""))

#generate plot
plot <- ggplot(avg_fields) + geom_point(aes(x=Actual, y = Predicted, 
                                           color=factor(Year), size=1),shape=1) + geom_abline(intercept=0, slope=1) +
  labs(x = "Actual yield (t/ha.)",y="Predicted yield (t/ha.)", title="Random Forest Actual Vs. Predicted Yield", color="Year") +
  scale_x_continuous(breaks = seq(0,60, by=5)) + scale_y_continuous(breaks = seq(0,60, by=5)) + coord_cartesian(ylim = c(27,55), xlim=c(27,55)) +
  theme(text = element_text(size = 18),legend.text=element_text(size=10)) +
  guides(size="none",colour = guide_legend(override.aes = list(size=5)))

ggsave(paste(dir,"/",dataname,"_actVpredPlot",run,".png",sep=""), width=7, height=5)

