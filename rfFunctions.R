#****************************************************************************
#Functions used in the RandomForest.R script. 
#See rfFunctions read.me for more information.
#****************************************************************************

#functions----------------------------------------------------------------------
calc_metrics <- function(actual,predicted){
  tryCatch(sd(predicted))
  mse <- mse(actual, predicted)
  rmse <- rmse(actual, predicted)
  rrmse <- rmse/mean(actual)
  r2 <- cor(actual,predicted)^2*100
  bias <- bias(actual, predicted)
  pbias <- bias/mean(actual)*100
  return (c(mse,rmse,rrmse,r2,bias,pbias))
}

min_mse <- function(fs_results){
  #Find best model based on MSE for add-on-in feature selection models
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
  #find the lowest test error while ensuring that the difference between train
  #and test error is not too large and is positive
  min_errors <- list(train_mse=min_train,train_vars=train_n,
                     test_mse=min_test,test_vars=test_n)
  a <- features[[min_errors$test_vars]]$metrics$train_mse
  b <- features[[min_errors$train_vars]]$metrics$test_mse
  
  diff_test <- min_errors$test_mse - a
  diff_train <- b - min_errors$train_mse
  
  if(diff_test < diff_train & diff_test > 0){
    return (fs_results[[test_n]])}else{return (fs_results[[train_n]])}
}

best_mod <- function(models, nodesize = F){
  #find best model based on RMSE
  if (!nodesize){
    mintrain <- 1000 #arbitrary initial value
    amod <- NULL
    mintest <- 1000 #arbitrary initial value
    bmod <- NULL
    
    btrain <- mod$train_metrics[2]
    bntree <- mod$rf$dots$ntree
    bmtry <- mod$rf$bestTune$mtry
    bpreds <- mod$prediction
    
    #find model with min train and test rmse, resp.
    for (mod in models){
      if(mod$test_metrics[2] < mintest){
        mintest <- mod$test_metrics[2]
        bmod <- mod
      }
      if(mod$train_metrics[2] < mintrain){
        mintrain <- mod$train_metrics[2]
        amod <- mod
      }
    }
    #find the lowest test error while ensuring that the difference between train
    #and test error is not too large and is positive
    diff_test <- mintest - bmod$train_metrics[2]
    diff_train <- amod$test_metrics[2] - mintrain
    
    if(diff_test < diff_train & diff_test > 0){
      bestmod <- bmod
    } else{bestmod <- amod}
    btest <-  bestmod$test_metrics[2]
    btrain <- bestmod$train_metrics[2]
    bntree <- bestmod$rf$dots$ntree
    bmtry <- bestmod$rf$bestTune$mtry
    bpreds <- bestmod$prediction
    return (list(ntree=bntree,mtry=bmtry,train_pred=bestmod$rf$finalModel$predicted,
                 predictions=bpreds,
                 test_metrics=list(mse=bestmod$test_metrics[1],
                                   rmse=bestmod$test_metrics[2],
                                   rrmse=bestmod$test_metrics[3],
                                   r2=bestmod$test_metrics[4],
                                   bias=bestmod$test_metrics[5],
                                   pbias=bestmod$test_metrics[6])))
  } 
  else{ #find best model from nodesize gridsearch
    btest <- 1000
    ns <- names(models)
    i <- 1
    for (mod in models) {
      if(mod$test_metrics$rmse < btest){
        btest <- mod$test_metrics$rmse
        bmod <- mod
        nsize <- ns[i]
      }
      i <- i + 1
    }
    return(list(nodesize=as.numeric(nsize),model=bmod))
  }
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
  while (j <= v + 2) { #try additional 2 features than what CV gives
    set.seed(seed)
    trnx <- train_x[importance$vars[1:j]]
    tstx <- test_x[importance$vars[1:j]]
    if (j < 3){
      n <- 1
    } else { 
      mtry <- as.data.frame(tuneRF(trnx, train_y, stepFactor = 1.5))
      n <- mtry$mtry[mtry$OOBError == min(mtry$OOBError)]
    }
    
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
hype_tune <- function(train, test, y, cvmethod='cv',folds=10, reps=2, mtry_seq=(1:30), 
                      ntree_seq=seq(100,1000,by=100),nodesize=NULL,tune_search='grid'){
  #train - all training data (including train y)
  #test - test X values 
  #y - test y values
  starttime <- Sys.time()
  modellist <- list()
  
  control <- trainControl(method=cvmethod, 
                          number=folds, 
                          repeats=reps, 
                          search=tune_search,
                          savePredictions = 'final')
  
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
                           ntree=nt,
                           nodesize=nodesize)
    train_res <- calc_metrics(rf_gridsearch$pred$obs,rf_gridsearch$pred$pred)
    
    rf_pred <- predict(rf_gridsearch,test)
    pred_res <- calc_metrics(y,as.numeric(rf_pred))
    
    n <- as.numeric(rf_gridsearch$bestTune$mtry)
    print(paste("best mtry:",n, 
                "train_RMSE:",rf_gridsearch$results$RMSE[rf_gridsearch$results$mtry==n], 
                "test_RMSE:", pred_res[2],sep=" "))
    print(Sys.time() - starttime)
    info <- list(rf = rf_gridsearch, train_metric=train_res, prediction=rf_pred, test_metrics=pred_res)
    modellist[[key]] <- info
    
  }
  return (modellist)
}
