#organize preprocessing metrics into one table
#Actual and predicted yields for baseline, US,OS, SMOTER
#Training and testing errors for baseline, US,OS, SMOTER: RMSE, RRMSE,PBIAS

load("MLalone_slope-lidar2.RData")

dir <- "/preProcessPlots"
orgd <- "/MR_RForgTrainSets"


dirs <- c("/MR_RFnoCumClim-cv2","/MR_RFUndersampB-cv2", 
          "/MR_RFOversampB-cv2","/MR_RFSMOTER_B-cv2")
runs <- c("baseline","US","OS","SMOTER")
summ <- read.csv(paste(getwd(),dir,"/summary.csv",sep=""))
mdf <- data.frame()
fields <- labels[c(1,7,8,9,10,12,16,17,18)]
#organize data
first <- T
j <- 1
for(d in dirs){
  curd <- paste(getwd(),d,"/",sep="")
  files <- list.files(curd,recursive = T)
  mfiles <- files[grepl("metrics.csv",files)]
  if(first){
    mfiles <- mfiles[c(1,7,8,9,10,12,16,17,18)]
    first <- F
  }
  metrics <- lapply(paste(curd,mfiles,sep=""),read.csv)
  i <- 1
  for(m in metrics){
    temp <- cbind("Run"= runs[j], "Field"=labels[i],m[1,1:3],m[2,1:3],m[3,1:3],m[6,1:3])
    mdf <- rbind(mdf,temp)
    i <- i + 1
  }
  j <- j + 1
  # if (first){
  #   values <- targd
  #   first <- F
  # } else {
  #   values <- merge(x = values,y=targd,by="CFIDYr",all.x=F)
  # }
  # 
  # i <- 1
  # nt <- data.frame(values$CFIDYr)
  # while(i <= length(fields)){
  #   if (length(grep(paste(fields[i],"/metrics",sep=""),files))!=0){
  #     metrics <- read.csv(paste(curd,files[grep(paste(fields[i],"/metrics",sep=""),files)],sep=""))
  #     
  #     nt$RRMSE[match(labels[i],nt$values.CFIDYr)] <- metrics[3,3]*100
  #   } else{nt$RRMSE[match(labels[i],nt$values.CFIDYr)] <- NA}
  #   i <- i + 1
  # }
  # colnames(nt) <- c("CFIDYr",paste(d,"_RRMSE",sep=""))
  # values <- merge(values,nt)
}
colnames(values) <- c("CFIDYr","Year","Actual","baseline.Predicted","baseline.RRMSE",
                      "Y1","A1","US.Predicted","US.RRMSE",
                      "Y2","A2","OS.Predicted","OS.RRMSE",
                      "Y3","A3","SMOTE.Predicted","SMOTE.RRMSE")
values <- values[order(values$Year),] %>% select(1:5,8,9,12,13,16,17)
valuest <- as.data.frame(t(values))
write.csv(values,paste(getwd(),dir,"/summary.csv",sep=""),row.names = F)

values <- read.csv(paste(getwd(),dir,"/summary.csv",sep=""))
colnames(values) <- c("CFIDYr","Year","Actual","baseline.Predicted","baseline.RRMSE",
                      "US.Predicted","US.RRMSE",
                      "OS.Predicted","OS.RRMSE",
                      "SMOTER.Predicted","SMOTER.RRMSE")