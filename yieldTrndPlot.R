packages <- c("dplyr","readxl","ggplot2","tidyr","reshape2","RColorBrewer","stringr")
lapply(packages, library, character.only=TRUE)

dir <- "/MR_RFnoNormcv2_2"

#gridded data
data <- read.csv(paste(getwd(),dir,"/allFieldsnoNormcv2_2.csv",sep=""))
data$year <- as.numeric(str_sub(data$Label,-4,-1))

yld_trnd <- ggplot(data, aes(x=year,y=Actual, group=year,fill=year)) + geom_boxplot() +
  theme_bw() + labs(x="Year", y="Yield (t/ha.)") + 
  scale_fill_gradientn(colours=RColorBrewer::brewer.pal(6,"Paired"), guide="none") +
  scale_x_continuous(breaks = seq(2015,2021, by=1)) + 
  scale_y_continuous(breaks = seq(25,65, by=5)) 
yld_trnd

ggsave(paste(getwd(),dir,"/yieldBoxPlot.png",sep=""), width=7, height=5)

#by field
yld_trnd <- ggplot(data, aes(x=year,y=Actual, group=Label,fill=year)) + geom_boxplot() +
  theme_bw() + labs(x="Year", y="Yield (t/ha.)") + 
  scale_fill_gradientn(colours=RColorBrewer::brewer.pal(6,"Paired"), guide="none") +
  scale_x_continuous(breaks = seq(2015,2021, by=1)) + 
  scale_y_continuous(breaks = seq(25,65, by=5)) 
yld_trnd
ggsave(paste(getwd(),dir,"/yieldBoxPlot_byFld.png",sep=""), width=7, height=5)

#field averages
datafld <- read.csv(paste(getwd(),dir,"/ActualVsPredicted_avgnoNormcv2_2_ver2.csv",sep=""))

fldyld_trnd <- ggplot(datafld, aes(x=Year,y=Actual, group=Year,fill=Year)) + geom_boxplot() +
  theme_bw() + labs(x="Year", y="Yield (t/ha.)")
fldyld_trnd
