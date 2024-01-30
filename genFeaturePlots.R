packages <- c("dplyr","readxl","randomForest","caret","Metrics","ggplot2",
              "tidyr","reshape2","RColorBrewer")
lapply(packages, library, character.only=TRUE)

#dir <- "/MR_RFnoNormcv2_2"
dir <- "/MR_Hybridcv2_2"
  
rootPath=paste(getwd(),dir,sep="")
listFiles=list.files(rootPath,recursive=TRUE,full.names=TRUE)

featFiles=grepl("features",listFiles)

all_features <- c() #list of all feature occurences for all field-years

listFiles <- listFiles[featFiles]
#listFiles <- listFiles[c(1:5,7,8,10,11,13,14,16,17,19,21,23,25,27,28,29)]
listFiles <- listFiles[c(1:5,7,8,10,11:13,15,16,18,20,21:25)] #hybrid

for (f in listFiles){
  featList <- read.csv(f)
  all_features <- append(all_features,featList[,2])
}

all_features <- as.data.frame(all_features)

all_vars <- unique(all_features)
all_vars$first <-0
all_vars$second <- 0
all_vars$third <- 0
all_vars$below_third <- 0

for (f in listFiles){
  featList <- read.csv(f)
  i <- match(featList[1,2],all_vars$all_features)
  all_vars$first[i] <- all_vars$first[i] + 1
  j <- match(featList[2,2],all_vars$all_features)
  all_vars$second[j] <- all_vars$second[j] + 1
  k <- match(featList[3,2],all_vars$all_features)
  all_vars$third[k] <- all_vars$third[k] + 1
  m <- 4
  while (m <= nrow(featList)){
    i <- match(featList[m,2],all_vars$all_features)
    all_vars$below_third[i] <- all_vars$below_third[i] + 1
    m <- m + 1
  }
}

#plot in order of frequency-----------------------------------------------------
cc <- scales::seq_gradient_pal("red", "blue", "Lab")(seq(0,1,length.out=22))

fplot <- ggplot(all_features,aes(y=fct_infreq(all_features), fill=fct_infreq(all_features))) +  geom_bar() +
  labs(x="Frequency", y="Features") + scale_fill_manual(values=cc, guide="none") + theme_bw()

fplot
ggsave(paste(getwd(),dir,"/featurePlot.png",sep=""), width=7, height=5)


#plot stacked with top 3--------------------------------------------------------
total <- c()
lvls <- c("first","second","third","below third")
n <- nrow(all_vars)
for (x in all_vars$all_features){
  total <- append(total,rep(x,4))
}
total <- as.data.frame(total)
total$importance <- rep(lvls,22)

total$frequency <- 0
i <- 1
while (i <= nrow(all_vars)){
  j <- match(all_vars$all_features[i],total$total)
  total$frequency[j] <- all_vars$first[i]
  total$frequency[j+1] <- all_vars$second[i]
  total$frequency[j+2] <- all_vars$third[i]
  total$frequency[j+3] <- all_vars$below_third[i]
  i <- i + 1
}
colnames(total) <- c("all_features","importance","frequency")

all_vars$sum <- all_vars$first + all_vars$second + all_vars$third + all_vars$below_third
total <- merge(total, all_vars, by=c("all_features")) 
total <- total[order(total$sum, decreasing = T),]

all_vars <- all_vars[order(all_vars$sum, decreasing = T),]
level_order <- all_vars$all_features
importance_order <- c("first","second","third","below third")

fplot <- ggplot(total, aes(fill=importance, y=frequency, x=all_features)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(limits = level_order) +
  scale_fill_discrete(limits = importance_order) +
  coord_flip() +
  labs(x = "Features", y="Frequency",fill="Position")
fplot
ggsave(paste(getwd(),dir,"/featurePlot_Top3.png",sep=""), width=7, height=5)
