#Code to clean data. Final dataset saved as .Rdata file

lapply(c("readxl","dplyr"), library, character.only=TRUE)


#For Hybrid Data--------------------------
data <- data.frame(read.csv(paste(getwd(),"/MattRamsay_RU_20152021/mrDatawSTICS.csv",sep="")))
ss_data <- subset(data, data$sh_PER_SAN != 0) %>% select(2,8:33,35,42:47,50:55,59,60,63:70,73:80,83:90,93:169)

#For ML Alone-----------------------------
data <- data.frame(read.csv(paste(getwd(),"/MattRamsay_RU_20152021/mrData.csv",sep="")))
#select only datapoints that have %sand info (have soil samples)
ss_data <- subset(data, data$sh_PER_SAN != 0) %>% select(6:33,35,41:47,50:55,58,59,62:69,72:79,82:89,92:108)

  
# yld_dist_snd <- ggplot(ss_data, aes(x=yield_tha)) + geom_histogram(color="black", fill="white") + labs(x="Yield (t/ha)", title = "Distribution of Russet Yield Datapoints with Field Soil Information")
# yld_dist_snd
ss_data$PH1[which(ss_data$PH1 %in% c(0))] <- 6.0 #replace missing pH values with pH of 6 as per CANSIS data

#move field ID, year, and yield columns to front
ss_data <- ss_data %>% dplyr::select("plant_ye", everything())
ss_data <- ss_data %>% dplyr::select("CFIDYr", everything())
ss_data <- ss_data %>% dplyr::select("yield_tha", everything()) 

c <- length(ss_data)

fields <- unique(ss_data$CFIDYr) %>% sort() #get list of full field names
labels <- c("L-2015","OS-2017","I17-2016","L-2017","CN-2017","P-2017","BR-2017",
            "JM-2017","M1-2017","I13-2017","OCR-2019","IC-2019","I17-2019",
            "RB-2019","P-2019","M4-2017","OS-2019","L-2021","JM-2020",
            "GRS-2019") %>% sort() #general labels for fields

#Add additional point data features: lidar elevation, rainLatene
all_data <- read_excel(paste(getwd(),"/MattRamsay_alldata/agData-v02152023.5.xlsx",sep=""), guess_max = 5000)

new_feats <- all_data[c("CFIDYr","DEMslope1","lidar_elev","rainLatene","sh_PER_SAN")]
ss_data <- merge(ss_data,new_feats,by=c("DEMslope1","CFIDYr")) %>% distinct()

rm("all_data","data","new_feats")
