#For Hybrid Data--------------------------
data <- data.frame(read.csv(hybrid_data))

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
