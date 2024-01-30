#****************************************************************************
#Generate metrics and Actual Vs. Predicted plot for results of random forest overall.
#****************************************************************************

generate_plot <- function(datafile,run,method,label="Year"){
avg_fields <- read.csv(datafile)

metrics_oa <- calc_metrics(avg_fields$Actual, avg_fields$Predicted)

metrics_all <- data.frame("Overall Field-Years"=metrics_oa)
rownames(metrics_all) <- c("MSE","RMSE","RRMSE","R2","BIAS","PBIAS")
write.csv(metrics_all,paste(method,"_metrics",run,".csv",sep=""))

#generate plot
if(label=="CFIDYr"){
  plot <- ggplot(avg_fields) + 
    geom_point(aes(x=Actual, y = Predicted, 
                   color=factor(CFIDYr), size=1),shape=1) +
    geom_abline(intercept=0, slope=1) +
    labs(x = "Actual yield (t/ha.)",y="Predicted yield (t/ha.)", 
         title="Random Forest Actual Vs. Predicted Yield", color="CFIDYr") +
    scale_x_continuous(breaks = seq(0,60, by=5)) + 
    scale_y_continuous(breaks = seq(0,60, by=5)) + 
    coord_cartesian(ylim = c(27,55), xlim=c(27,55)) +
    theme(text = element_text(size = 18),legend.text=element_text(size=10)) +
    guides(size="none",colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(paste(method,"_actVpredCFID",run,".png",sep=""), width=7, height=5)
}
else{
  plot <- ggplot(avg_fields) + 
    geom_point(aes(x=Actual, y = Predicted,
                   color=factor(Year), size=1),shape=1) + 
    geom_abline(intercept=0, slope=1) +
    labs(x = "Actual yield (t/ha.)",y="Predicted yield (t/ha.)", 
         title="Random Forest Actual Vs. Predicted Yield", color="Year") +
    scale_x_continuous(breaks = seq(0,60, by=5)) + 
    scale_y_continuous(breaks = seq(0,60, by=5)) + 
    coord_cartesian(ylim = c(27,55), xlim=c(27,55)) +
    theme(text = element_text(size = 18),legend.text=element_text(size=10)) +
    guides(size="none",colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(paste(method,"_actVpredYear",run,".png",sep=""), width=7, height=5)
}

}