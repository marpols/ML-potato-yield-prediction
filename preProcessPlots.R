library(knitr)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggh4x)

load("MLalone_slope-lidar2.RData")

dir <- "/preProcessPlots"
dir.create(paste(getwd(),dir,sep=""))
orgd <- "/MR_RForgTrainSets"


dirs <- c("/MR_RFnoCumClim-cv2","/MR_RFUndersampB-cv2", 
          "/MR_RFOversampB-cv2","/MR_RFSMOTER_B-cv2")

#organize data
first <- T
for(d in dirs){
  curd <- paste(getwd(),d,"/",sep="")
  files <- list.files(curd,recursive = T)
  targd <- read.csv(paste(curd,files[1],sep=""))

  if (first){
    values <- targd
    first <- F
  } else {
    values <- merge(x = values,y=targd,by="CFIDYr",all.x=F)
  }
  
  i <- 1
  nt <- data.frame(values$CFIDYr)
  while(i <= length(fields)){
    if (length(grep(paste(fields[i],"/metrics",sep=""),files))!=0){
      metrics <- read.csv(paste(curd,files[grep(paste(fields[i],"/metrics",sep=""),files)],sep=""))

      nt$RRMSE[match(labels[i],nt$values.CFIDYr)] <- metrics[3,3]*100
    } else{nt$RRMSE[match(labels[i],nt$values.CFIDYr)] <- NA}
    i <- i + 1
  }
  colnames(nt) <- c("CFIDYr",paste(d,"_RRMSE",sep=""))
  values <- merge(values,nt)
}
colnames(values) <- c("CFIDYr","Year","Actual","Baseline.Predicted","Baseline.RRMSE",
                      "Y1","A1","US.Predicted","US.RRMSE",
                      "Y2","A2","OS.Predicted","OS.RRMSE",
                      "Y3","A3","SMOTE.Predicted","SMOTE.RRMSE")
values <- values[order(values$Year),] %>% select(1:5,8,9,12,13,16,17)
valuest <- as.data.frame(t(values))
write.csv(values,paste(getwd(),dir,"/summary.csv",sep=""),row.names = F)

values <- read.csv(paste(getwd(),dir,"/summary.csv",sep=""))
colnames(values) <- c("CFIDYr","Year","Actual","Baseline.Predicted","Baseline.RRMSE",
                      "US.Predicted","US.RRMSE",
                      "OS.Predicted","OS.RRMSE",
                      "SMOTER.Predicted","SMOTER.RRMSE")
valuest <- as.data.frame(t(values))


names <- c(rep("Field 1a",9),rep("Field 2",9),rep("Field 3a",9),rep("Field 4a",9),
           rep("Field 5",9),rep("Field 6",9),rep("Field 4b",9),rep("Field 3b",9), rep("Field 1b",9))
years <- sort(rep(values$Year,9))

#Fix sort so that it is by factor, not alphabetical
plotvals <- data.frame(names,years,c(rep(row.names(valuest)[3:11],9)),
                       c(valuest[3:11,1],valuest[3:11,2],valuest[3:11,3],
                         valuest[3:11,4],valuest[3:11,5],valuest[3:11,6],
                         valuest[3:11,7],valuest[3:11,8],valuest[3:11,9]))
colnames(plotvals) <- c("Field","Year","Category","Value")

cats <- c("Baseline","US","OS","SMOTER")

errors <- plotvals[grepl("RRMSE",plotvals$Category),]
yields <- anti_join(plotvals,errors)

for(c in cats){
  errors[grepl(c,errors$Category),]$Category <- c
}
errors$Category <- factor(errors$Category, levels = unique(errors$Category),ordered = TRUE)
errors$Field <- factor(errors$Field, levels = c("Field 1a","Field 2","Field 3a","Field 4a","Field 5",
                                                "Field 6","Field 4b","Field 3b", "Field 1b"))
for(c in cats){
  yields[grepl(c,yields$Category),]$Category <- c
}
yields$Category <- factor(yields$Category, levels = unique(yields$Category),ordered = TRUE)
yields$Field <- factor(yields$Field, levels = c("Field 1a","Field 2","Field 3a","Field 4a","Field 5",
                                               "Field 6","Field 4b","Field 3b", "Field 1b"))

yields$Value <- as.numeric(yields$Value)
errors$Value <- as.numeric(errors$Value)



#Actual yield Vs. no pre-process Vs. under sampling Vs. Over sampling Vs. SMOTE predicted yield
#bar plot-----------------------------------------------------------------------
yields$Category <- forcats::fct_inorder(yields$Category)

bar_plot <- ggplot(yields, aes(x = Field, y = Value, fill = Category)) +
  geom_bar(stat = "identity",position = "dodge", width = 1.0) +
  labs(title = "Average Whole Field Yields: Actual and Predicted",
       x = " ",
       y = "Yield (t./ha.)",
       fill = " ") +
  theme_minimal() + 
  theme(legend.position = "bottom", axis.text.x=element_blank(), 
        legend.text=element_text(size=8), legend.key.size = unit(0.5,"cm")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")) +
  facet_wrap(~ Field + Year, scales="free_x",nrow=1, strip.position = "bottom") +
  scale_y_continuous(breaks=seq(0,60,by=5)) + coord_cartesian(ylim=c(0,60))

bar_plot
ggsave(paste(getwd(),dir,"/yieldsAvP4.jpg",sep=""), width=15, height=10, unit="cm")

#no pre-process Vs. under sampling Vs. Over sampling Vs. SMOTE predicted RRMSE
#bar plot ----------------------------------------------------------------------

errors$Category <- forcats::fct_inorder(errors$Category)
bar_plot <- ggplot(errors, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Errors of Field Predictions",
       x = " ",
       y = "RRMSE(%)",
       fill=" ") +
  theme_minimal() + 
  theme(legend.position = "bottom",axis.text.x=element_blank(),
        legend.key.size = unit(0.5,"cm")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73","#CC79A7")) +
  facet_wrap(~ Field + Year, scales="free_x",nrow=1, strip.position = "bottom") +
  scale_y_continuous(breaks=c(5,10,15,20,25,30,35,40))
bar_plot
ggsave(paste(getwd(),dir,"/rrmse_bp4.jpg",sep=""),width=15, height=10, unit="cm")


#Scatter Plot to show synthetic data points created by SMOTE for 2020-----------
orgdata <- read.csv(paste(getwd(),"/MR_RForgTrainSets/JimMac25-2020/data_orgTrainSets.csv",sep=""))
newdata <- read.csv(paste(getwd(),dirs[4],"/JimMac25-2020/newtrain_SMOTER_B-cv2.csv",sep=""))
synthetic <- anti_join(newdata,orgdata)
undersamp <- anti_join(orgdata,newdata)
opts <- semi_join(orgdata,newdata) #unchanged datapoints

synthetic$SMOTER <- T
newdata <- as.data.frame(merge(x=newdata,y=synthetic,all.x=T))
newdata$Yield_Category <- NA
newdata[newdata$yield_tha <= 38,]$Yield_Category <- "low" #low
newdata[newdata$yield_tha > 38 & newdata$yield_tha <= 49,]$Yield_Category <- "middle" #middle
newdata[newdata$yield_tha > 49,]$Yield_Category <- "high" #high

newdata.n<- as.data.frame(scale(newdata[,2:76]))

# newdata <- as.data.frame(scale(merge(x=newdata,y=synthetic,all.x=T)))
newdata.n[!is.nan(newdata.n$SMOTER),]$SMOTER <- "Original"
newdata.n[newdata.n$SMOTER == "NaN",]$SMOTER <- "SMOTER"

newdata[is.na(newdata$SMOTER),]$SMOTER <- "Original"
newdata[newdata$SMOTER==T,]$SMOTER <- "SMOTER"

newdata.n$Yield_Category <- newdata$Yield_Category



newdata.n$Yield_Category <- factor(newdata.n$Yield_Category, levels = c("low","middle","high"))
newdata$Yield_Category <- factor(newdata$Yield_Category, levels = c("low","middle","high"))

write.csv(newdata.n,paste(getwd(),dir,"/scaledData.csv",sep=""))


newdata.n$SMOTER <- factor(newdata.n$SMOTER, levels = c("Original","SMOTER"))

#scatterplot with colour for synthetic and shape for low, high yield

shuf.data <- sample(newdata)

ggplot(newdata.n, 
       aes(x = June.precip, y = June.GDD, color=SMOTER,shape=Yield_Category)) +
  geom_point(size = 4) +
  labs(title = "Scatter Plot of Differences",
       x = "June Precipitation",
       y = "June GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_area(breaks = c(1,2)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))
ggsave(paste(getwd(),dir,"/SMOTERplot1-forLgnd.png",sep=""), width=15, height=10, unit="cm")

ggplot(newdata.n, aes(x = July.precip, y = July.GDD, color=SMOTER,shape=Yield_Category, size=SMOTER)) +
  geom_point(stroke=1) +
  labs(title = "Scatter Plot of Differences",
       x = "July Precipitation",
       y = "July GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_manual(values= c(2,5)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))
ggsave(paste(getwd(),dir,"/SMOTERplot1-JM2020_JulyClim2_v2.png",sep=""), width=15, height=10, unit="cm")

# ggplot(newdata, aes(x = DEMslope1 , y = yield_tha, color=SMOTE)) +
#   geom_point() +
#   labs(title = "Scatter Plot of Differences",
#        x = "PC3",
#        y = "PC4",
#        color = " ") +
#   theme_minimal()
# ggsave(paste(getwd(),dir,"/SMOTERplot2-JM2020.jpg",sep=""))

#FREQUECY BARPLOT OF POINTS TO GO ABOVE SCATTER PLOT
x <- newdata.n[newdata.n$SMOTER=="SMOTER",] |> group_by(July.precip,July.GDD)
c <- count(x)
ggplot(c, aes(x=July.precip,y=n, fill="#00BFC4")) +
  geom_bar(stat = "identity") +
  labs(title = "Scatter Plot of Differences",
       x = "June Precipitation",
       y = "Number of points",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_fill_manual(values=c("#00BFC4"))
ggsave(paste(getwd(),dir,"/SMOTERplot1-JM2020_JulyClim2_freq.png",sep=""), width=15, height=5, unit="cm")

## 3 dimensions
colors <- setNames(rainbow(length(unique(newdata$SMOTE))), unique(newdata$SMOTE))

splot <- scatterplot3d(newdata$DEMslope1, newdata$sh_PER_OM, newdata$yield_tha,
              main = "3D Scatter Plot", color = colors[newdata$SMOTE], pch = 20, angle=55,
              xlab="PC1",ylab="PC2",zlab="PC3",)

#--Histogram of under/oversampled data from SMOTER-------------

synthetic$Category <- "Synthetic"
undersamp$Category <- "Removed"
opts$Category <- "Original"
s <- data.frame("Yield" = synthetic$yield_tha,"Category"= synthetic$Category) #synthetic points
u <- data.frame("Yield" = undersamp$yield_tha,"Category"= undersamp$Category) #undersampled points
o <- data.frame("Yield" = opts$yield_tha,"Category"= opts$Category) #original points

catd <- rbind(o,u,s) #new data frame
catd$Yield <- as.double(catd$Yield)

catd$Category <- factor(catd$Category, levels = c("Removed","Synthetic","Original"))


#all categories w/ legend
ggplot(catd, aes(x=Yield, fill=Category)) + 
  geom_histogram(position="stack",bins=26) +
  labs(title = "Training Dataset Preprocessed with SMOTER",
       x = "Yield (t./ha.)",
       y = "Count",
       fill = " ") +
  theme_minimal() + 
  guides() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("grey","#00BFC4", "#F8766D"))
ggsave(paste(getwd(),dir,"/SMOTERtrain-JM20202forlengend.png",sep=""), height=10, width=15, unit="cm")

#synthetic and original, no legend
ggplot(catd[catd$Category != "Removed",], aes(x=Yield, fill=Category)) + 
  geom_histogram(position="stack",bins=26) +
  labs(title = "Training Dataset Preprocessed with SMOTER",
       x = "Yield (t./ha.)",
       y = "Count",
       fill = " ") +
  theme_minimal() + 
  guides(fill = "none") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("#00BFC4", "#F8766D"))
ggsave(paste(getwd(),dir,"/SMOTERtrain-JM20202_11.png",sep=""), height=10, width=15, unit="cm")

#original only, white fill, black outline
ggplot(catd[catd$Category=="Original",], aes(x=Yield, fill=Category)) + 
  geom_histogram(position="stack",bins=26,color="black",size=1) +
  geom_histogram(bins=26) +
  labs(title = "Training Dataset Preprocessed with SMOTER",
       x = "Yield (t./ha.)",
       y = "Count",
       fill = " ") +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("white","#00BFC4", "#F8766D"))
ggsave(paste(getwd(),dir,"/forlegend2.png",sep=""), height=10, width=15, unit="cm")

#original, outlines
ggplot(orgdata,aes(x=yield_tha)) + 
  geom_histogram(bins=26, color="black", size=1, fill="white",alpha=0) +
  labs(title = "Training Dataset Preprocessed with SMOTER",
       x = "Yield (t./ha.)",
       y = "Count",
       fill = " ") +
  theme_minimal() + 
  guides(fill = "none") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste(getwd(),dir,"/SMOTERtrain-JM20202_11.png",sep=""), height=10, width=15, unit="cm")

ggplot(catd[catd$Category != "synthetic",], aes(x=Yield, fill=Category)) + 
  geom_histogram(position="stack",bins=26, color="black", size=1.5) +
  geom_histogram(bins=26, fill="white") +
  labs(title = "Training Dataset Preprocessed with SMOTER",
       x = "Yield (t./ha.)",
       y = "Count",
       fill = " ") +
  theme_minimal() + 
  guides(fill = "none") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste(getwd(),dir,"/SMOTERtrain-JM20202_8.png",sep=""), height=10, width=15, unit="cm")

ggplot(orgdata,aes(x=yield_tha)) +
  geom_histogram(fill="#FF9999") +
  labs(title = "Original Training Dataset",
       x = "Yield (t./ha.)",
       y = "Count") +
  theme_minimal()
ggsave(paste(getwd(),dir,"/Orgtrain-JM2020.jpg",sep=""))



#boxplots
x1 <- read.csv(paste(getwd(),"/MR_RForgTrainSets/JimMac25-2020/data_orgTrainSets.csv",sep="")) |> select(2)
x2 <- read.csv(paste(getwd(),dirs[4],"/JimMac25-2020/newtrain_SMOTER_B-cv2.csv",sep=""))|> select(2)
n <- c("yield","category")

x1$Category <- "Baseline"
x2$category <- "SMOTER"
colnames(x1) <- n
colnames(x2) <- n

dfbox <- rbind(x1,x2)
colnames(dfbox) <- n
 
ggplot(dfbox,aes(x=yield, y=category)) +
  geom_boxplot() +
  labs(title = " ",y = " ", x = "Yield (t./ha.)") +
  theme_minimal()
ggsave(paste(getwd(),dir,"/bx-JM2020.jpg",sep=""),height=10, width=15, unit="cm")



#----------------------- Soil and Management Table------------------------------
soil_info <- read.csv("C:/Users/PolsinelliM/OneDrive - AGR-AGR/Documents/Potato Yield Prediction/ML/MattRamsay_RU_20152021/soil_fieldInfo.csv")
soil_info <- subset(soil_info,select=-1) %>%
  mutate_if(is.numeric, ~round(., 3))

colnames(soil_info) <- c("","Year","% Sand","% Silt", "% Clay", "% OM", "FC", 
                          "WP", "BD","PAW", "Sowing", "Harvest", "Yield")
soil_info$Year <- factor(soil_info$Year, 
                         levels = c("2015","2017","2019","2020","2021"),
                         ordered = TRUE)


table <- kbl(soil_info[c(1,7,8,9,10,12,16,17,18),],row.names = FALSE) %>% kable_paper("striped", full_width = F) %>% 
  add_header_above(c(" "=2, "Soil Texture" = 4, "Field Attributes" = 4, " "=3))



#----------------------- Climate Plot (GDD and Precip) -------------------------
cdf <- read.csv("C:/Users/PolsinelliM/OneDrive - AGR-AGR/Documents/Potato Yield Prediction/ML/Climate Data/cum30yr_climateData.csv")

summ <- summarise(cdf,.by=Year,"GDD"=sum(GDD),"Precipitation"=sum(Precip.),"Cat"="Yearly")
tya_precip <- mean(summ[9:39,3])
tya_gdd <- mean(summ[9:39,2])
df <- data.frame("Year"=summ$Year,"GDD"=tya_gdd,"Precipitation"=tya_precip,"Cat"="30 Year Avg.")
summ$tyaGDD <- tya_gdd
summ$tyaPrecip <- tya_precip
summ2 <- rbind(summ,df)


ggplot(summ2[summ2$Year %in% c(2015,2017,2019,2020,2021),], aes(x=Year, y=GDD, color=Cat)) + 
  geom_line(aes(linetype = Cat, size = Cat)) + 
  geom_point(data = summ2[summ2$Cat=="Yearly"&summ2$Year %in% c(2015,2017,2019,2020,2021),], size=2) + 
  theme_minimal() + guides(color="none", size="none", linetype="none") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_size_manual(values = c(1.5,1)) + 
  labs(y="Cu GDD (Base 5)", color=" ", size =" ", linetype = " ") +
  scale_y_continuous(breaks = seq(1650,2600, by=50)) +
  scale_x_continuous(breaks = c(2015,2017,2019,2020,2021))
ggsave(paste(getwd(),dir,"/yearlyGDD3.jpg",sep=""),width=2000, height=625, unit="px")

ggplot(summ2[summ2$Year %in% c(2015,2017,2019,2020,2021),], aes(x=Year, y=Precipitation, color=Cat)) + 
  geom_line(aes(linetype = Cat, size = Cat)) + 
  geom_point(data = summ2[summ2$Cat=="Yearly"&summ2$Year %in% c(2015,2017,2019,2020,2021),], size=2) + 
  theme_minimal() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_size_manual(values = c(1.5,1)) + 
  labs(y="Precipiation (mm)", color=" ", size =" ", linetype = " ") +
  scale_y_continuous(breaks = seq(650,1150, by=50)) +
  scale_x_continuous(breaks = c(2015,2017,2019,2020,2021))
ggsave(paste(getwd(),dir,"/yearlyPrecip3.jpg",sep=""),width=2000, height=875, unit="px")


#------Actual Vs. Predicted Yield for 9 fields----------------------------------
#Scatter plot to visualize overall performance
svals <- data.frame("Field" = c("Field 1a","Field 2","Field 3a","Field 4a","Field 5",
                                "Field 6","Field 4b","Field 3b", "Field 1b"),
                    "Actual" = yields[order(yields$Category),][1:9,c(2,4)],
                    "Baseline" = yields[order(yields$Category),][10:18,4],
                    "US" = yields[order(yields$Category),][19:27,4],
                    "OS" = yields[order(yields$Category),][28:36,4],
                    "SMOTER" = yields[order(yields$Category),][37:45,4])
svals$Actual.Year <- as.character(svals$Actual.Year)
svals$Actual.Value <- as.numeric(svals$Actual.Value)
svals$Baseline <- as.numeric(svals$Baseline)
svals$US <- as.numeric(svals$US)
svals$OS <- as.numeric(svals$OS)
svals$SMOTER <- as.numeric(svals$SMOTER)


plot <- ggplot(svals, aes(x=Actual.Value, y = Baseline, color=Actual.Year, shape=Actual.Year)) + 
  geom_point(size=5, stroke = 1.2) + theme_minimal() +
  geom_abline(intercept=0, slope=1) +
  labs(x = "Actual yield (t./ha.)",y="Predicted yield (t./ha.)", 
       title="Baseline", color=" ", shape = " ") +
  scale_x_continuous(breaks = seq(0,60, by=5)) + 
  scale_y_continuous(breaks = seq(0,60, by=5)) + 
  coord_cartesian(ylim = c(27,55), xlim=c(27,55)) +
  theme(text = element_text(size = 18),legend.text=element_text(size=10)) +
  guides(size="none",colour = guide_legend(override.aes = list(size=5))) +
  scale_shape_manual(values = c(0,1,2,3,4))

plot
ggsave(paste(getwd(),dir,"/ActVPred_baseline3.jpg",sep=""), width=7, height=5)
df <- data.frame("RMSE"= Metrics::rmse(svals$Actual.Value,svals$OS),
                 "RRMSE" = Metrics::rmse(svals$Actual.Value,svals$OS)/mean(svals$Actual.Value)*100,
                "MAE"= Metrics::mae(svals$Actual.Value,svals$OS),
                "PBIAS"= Metrics::percent_bias(svals$Actual.Value,svals$OS))
write.csv(df,paste(getwd(),dir,"/OSmetrics.csv",sep=""))

#-----train test RRMSE bar graphs-----------------------------------------------
ttdf <- read_excel(paste(getwd(),dir,"/metricsSummary.xlsx",sep=""), sheet="Sheet3") #train/test set errors
ttdf$Set <- factor(ttdf$Set, levels=c("Train","Test"))
ttdf$Category <- factor(ttdf$Category, levels=c("Baseline","US","OS","SMOTER"))


#X-axis = Field, Y-axis = Category
ggplot(ttdf, aes(x = Set, y = RRMSE, fill = Set)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  labs(title = "Errors of Field Predictions",
       x = " ",
       y = "RRMSE(%)",
       fill=" ") +
  theme_minimal() + 
  theme(legend.position = "right",axis.text.x=element_blank(),
        legend.key.size = unit(0.5,"cm"),
        strip.placement = "outside") +
  scale_fill_manual(values=c("#56B4E9", "#009E73")) +
  facet_grid(Category ~ Field,switch="x", scale="free") +
  force_panelsizes(rows=1, cols=c(0.3,0.3,0.3,0.3))

ggsave(paste(getwd(),dir,"/trainTest_RRMSEfr.jpg",sep=""), width=23, height=10, unit="cm")

#X-axis = Category, Y-axis = Field
ggplot(ttdf, aes(x = Set, y = RRMSE, fill = Set)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  labs(title = "Errors of Field Predictions",
       x = " ",
       y = "RRMSE(%)",
       fill=" ") +
  theme_minimal() + 
  theme(legend.position = "right",axis.text.x=element_blank(),
        legend.key.size = unit(0.5,"cm"),
        strip.placement = "outside") +
  scale_fill_manual(values=c("#56B4E9", "#009E73")) +
  facet_grid(Field ~ Category,switch="x", scale="free") +
  force_panelsizes(rows=1, cols=c(0.3,0.3,0.3,0.3))

ggsave(paste(getwd(),dir,"/trainTest_RRMSEfr2.jpg",sep=""), width=10, height=20, unit="cm")
ggsave(paste(getwd(),dir,"/trainTest_RRMSEfr.jpg",sep=""), width=10, height=20, unit="cm")

#Colour by Category
ggplot(ttdf, aes(x = Set, y = RRMSE, fill = Category, alpha = Set)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  labs(title = "Errors of Field Predictions",
       x = " ",
       y = "RRMSE(%)",
       fill=" ") +
  theme_minimal() + 
  theme(legend.position = "right",
        legend.key.size = unit(0.5,"cm"),
        strip.placement = "outside") +
  guides(alpha = "none") +
  facet_grid(Field ~ Category,switch="x", scale="free") +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73","#CC79A7"))+
  scale_alpha_manual(values=c(0.5,1.5)) +
  force_panelsizes(rows=1, cols=c(0.3,0.3,0.3,0.3))

ggsave(paste(getwd(),dir,"/trainTest_RRMSEfr3.jpg",sep=""), width=15, height=23, unit="cm")


#----Analysing Other fields-----------------------------------------------------
#"Leitas-2015"
#'Leitas-2021"

fy <- "JimMac25-2017"

orgdata <- read.csv(paste(getwd(),"/MR_RForgTrainSets/",fy,"/data_orgTrainSets.csv",sep=""))
newdata <- read.csv(paste(getwd(),dirs[4],"/",fy,"/newtrain_SMOTER_B-cv2.csv",sep=""))
synthetic <- anti_join(newdata,orgdata)
undersamp <- anti_join(orgdata,newdata)
opts <- semi_join(orgdata,newdata) #unchanged datapoints

f <- ss_data[ss_data$CFIDYr == 'Leitas-2021',]
testset <- ss_data[ss_data$CFIDYr == fy,]
synthetic$SMOTER <- T

newdata <- as.data.frame(merge(x=newdata,y=synthetic,all.x=T))
newdata$Yield_Category <- NA
newdata[newdata$yield_tha <= 38,]$Yield_Category <- "low" #low
newdata[newdata$yield_tha > 38 & newdata$yield_tha <= 49,]$Yield_Category <- "middle" #middle
newdata[newdata$yield_tha > 49,]$Yield_Category <- "high" #high



newdata.n<- as.data.frame(scale(newdata[,2:76]))

# newdata <- as.data.frame(scale(merge(x=newdata,y=synthetic,all.x=T)))
newdata.n[!is.nan(newdata.n$SMOTER),]$SMOTER <- "Original"
newdata.n[newdata.n$SMOTER == "NaN",]$SMOTER <- "SMOTER"

newdata[is.na(newdata$SMOTER),]$SMOTER <- "Original"
newdata[newdata$SMOTER==T,]$SMOTER <- "SMOTER"

newdata.n$Yield_Category <- newdata$Yield_Category

newdata.n$Yield_Category <- factor(newdata.n$Yield_Category, levels = c("low","middle","high"))
newdata$Yield_Category <- factor(newdata$Yield_Category, levels = c("low","middle","high"))

#write.csv(newdata.n,paste(getwd(),dir,"/scaledData.csv",sep=""))

newdata.n$SMOTER <- factor(newdata.n$SMOTER, levels = c("Original","SMOTER"))

#scatterplot with colour for synthetic and shape for low, high yield

#June-------------------
ggplot(newdata, aes(x = June.precip, y = June.GDD, color=SMOTER,shape=Yield_Category, size=SMOTER)) +
  geom_point(stroke=1) +
  geom_abline(intercept=mean(f$June.GDD), slope=0) +
  geom_vline(xintercept=mean(f$June.precip)) +
  labs(title = "Scatter Plot of Differences",
       x = "June Precipitation",
       y = "June GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_manual(values= c(2,5)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))

ggplot(newdata, aes(x = June.precip, y = June.GDD, color=SMOTER,shape=Yield_Category, size=SMOTER)) +
  geom_point(stroke=1) +
  geom_abline(intercept=mean(testset$June.GDD), slope=0) +
  geom_vline(xintercept=mean(testset$June.precip)) +
  labs(title = "Scatter Plot of Differences",
       x = "June Precipitation",
       y = "June GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_manual(values= c(2,5)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))

ggsave(paste(getwd(),dir,"/SMOTERField1a-June.png",sep=""), width=15, height=10, unit="cm")
#frequency barplot
x <- newdata.n[newdata.n$SMOTER=="SMOTER",] |> group_by(June.precip,June.GDD)
c <- count(x)
ggplot(c, aes(x=June.precip,y=n, fill="#00BFC4")) +
  geom_bar(stat = "identity") +
  labs(title = "Scatter Plot of Differences",
       x = "June Precipitation",
       y = "Number of points",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_fill_manual(values=c("#00BFC4"))
ggsave(paste(getwd(),dir,"/SMOTERField1a_Junefreq.png",sep=""), width=15, height=5, unit="cm")

#July-----------------
ggplot(newdata, aes(x = July.precip, y = July.GDD, color=SMOTER,shape=Yield_Category, size=SMOTER)) +
  geom_point(stroke=1) +
  geom_abline(intercept=mean(f$July.GDD), slope=0) +
  geom_vline(xintercept=mean(f$July.precip)) +
  labs(title = "Scatter Plot of Differences",
       x = "July Precipitation",
       y = "July GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_manual(values= c(2,5)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))

ggplot(newdata, aes(x = July.precip, y = July.GDD, color=SMOTER,shape=Yield_Category, size=SMOTER)) +
  geom_point(stroke=1) +
  geom_abline(intercept=mean(testset$July.GDD), slope=0) +
  geom_vline(xintercept=mean(testset$July.precip)) +
  labs(title = "Scatter Plot of Differences",
       x = "July Precipitation",
       y = "July GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_manual(values= c(2,5)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))

ggsave(paste(getwd(),dir,"/SMOTERField1a-July.png",sep=""), width=15, height=10, unit="cm")
#frequency barplot
x <- newdata.n[newdata.n$SMOTER=="SMOTER",] |> group_by(July.precip,July.GDD)
c <- count(x)
ggplot(c, aes(x=July.precip,y=n, fill="#00BFC4")) +
  geom_bar(stat = "identity") +
  labs(title = "Scatter Plot of Differences",
       x = "July Precipitation",
       y = "Number of points",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_fill_manual(values=c("#00BFC4"))
ggsave(paste(getwd(),dir,"/SMOTERField1_Julyfreq.png",sep=""), width=15, height=5, unit="cm")

#August-----------------
ggplot(newdata, aes(x = Aug.precip, y = Aug.GDD, color=SMOTER,shape=Yield_Category, size=SMOTER)) +
  geom_point(stroke=1) +
  geom_abline(intercept=mean(f$Aug.GDD), slope=0) +
  geom_vline(xintercept=mean(f$Aug.precip)) +
  labs(title = "Scatter Plot of Differences",
       x = "August Precipitation",
       y = "August GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_manual(values= c(2,5)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))

ggplot(newdata, aes(x = Aug.precip, y = Aug.GDD, color=SMOTER,shape=Yield_Category, size=SMOTER)) +
  geom_point(stroke=1) +
  geom_abline(intercept=mean(testset$Aug.GDD), slope=0) +
  geom_vline(xintercept=mean(testset$Aug.precip)) +
  labs(title = "Scatter Plot of Differences",
       x = "August Precipitation",
       y = "August GDD",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_shape_manual(values = c(1,0,2)) +
  scale_size_manual(values= c(2,5)) +
  scale_fill_manual(values=c("#F8766D","#00BFC4"))

ggsave(paste(getwd(),dir,"/SMOTERField1a-August.png",sep=""), width=15, height=10, unit="cm")
#frequency barplot
x <- newdata.n[newdata.n$SMOTER=="SMOTER",] |> group_by(Aug.precip,Aug.GDD)
c <- count(x)
ggplot(c, aes(x=Aug.precip,y=n, fill="#00BFC4")) +
  geom_bar(stat = "identity") +
  labs(title = "Scatter Plot of Differences",
       x = "August Precipitation",
       y = "Number of points",
       color = " ",
       shape = "Yield") +
  theme_minimal()  +
  scale_fill_manual(values=c("#00BFC4"))
ggsave(paste(getwd(),dir,"/SMOTERField1a_Augfreq.png",sep=""), width=15, height=5, unit="cm")

#Actual V Predicted for gridded data
predSM <- read.csv(paste(getwd(),dirs[4],"/IansCornerE1-2019/preds_SMOTER_B-cv2.csv",sep=""))
predOS <- read.csv(paste(getwd(),dirs[3],"/IansCornerE1-2019/preds_OversampB-cv2.csv",sep=""))
predB <- read.csv(paste(getwd(),dirs[1],"/IansCornerE1-2019/preds_noCumClim-cv2.csv",sep=""))

plot <- ggplot(predB, aes(x=Actual, y = Predicted,)) + 
  geom_point(size=5, stroke = 1.2, shape=1) + theme_minimal() +
  geom_abline(intercept=0, slope=1) +
  labs(x = "Actual yield (t./ha.)",y="Predicted yield (t./ha.)", 
       title="Baseline") +
  scale_x_continuous(breaks = seq(0,60, by=5)) + 
  scale_y_continuous(breaks = seq(0,60, by=5)) + 
  coord_cartesian(ylim = c(32,42), xlim=c(32,42)) +
  theme(text = element_text(size = 18),legend.text=element_text(size=10))
plot

plot <- ggplot(predOS, aes(x=Actual, y = Predicted,)) + 
  geom_point(size=5, stroke = 1.2, shape=1) + theme_minimal() +
  geom_abline(intercept=0, slope=1) +
  labs(x = "Actual yield (t./ha.)",y="Predicted yield (t./ha.)", 
       title="OS") +
  scale_x_continuous(breaks = seq(0,60, by=5)) + 
  scale_y_continuous(breaks = seq(0,60, by=5)) + 
  coord_cartesian(ylim = c(32,42), xlim=c(32,42)) +
  theme(text = element_text(size = 18),legend.text=element_text(size=10))
plot
