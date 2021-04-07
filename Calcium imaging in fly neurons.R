# ANALYSIS OF TDTOMATO-GCAMP6 DATA

library(tidyverse)
library(readr)
library(MESS)
library(broom)

# Data file to be used should be a CSV formatted as indicated here.
# Column 1 should contain the numbers 1, 2, 3.. that mark the events such as removal of Ca2+ or GPN addition. This column should NOT have a label.
# Column 2 should contain time in seconds. This column can have any label.
# Column 3 onward should first contain the GFP intensities labeled GFP1, GFP2, and so on. The last of these columns should contain the background values for GFP.
# After the GFP intensities, the columns should contain the mCherry intensities, labeled mCherry1, mCherry2, and so on. The last of these columns should contain the background values for mCherry.
# Ratios from the beginning to the first addition event (1 in column 1) are considered the baseline and will be set to 1.

ratio_file <- read_csv("~/ratio.csv") # The name of the input file can be modified as needed

# The following block will separate the GFP and mCherry columns and subtract backgrounds for each.

col_names <- as.data.frame(colnames(ratio_file))
rem_cols <- which(grepl("mCherry",col_names$`colnames(ratio_file)`)==T)

ratio_file_GFP <- as.data.frame(ratio_file[,-c(rem_cols)])
for (j in 3:(ncol(ratio_file_GFP)-1)){
  col1 <- as.matrix(ratio_file_GFP[,c(j)])
  col_back <- as.matrix(ratio_file_GFP[,c(ncol(ratio_file_GFP))])
  ratio_file_GFP[,c(j)]<-col1-col_back
}
ratio_file_GFP[,c(ncol(ratio_file_GFP))]<-NULL

ratio_file_mCherry <- as.data.frame(ratio_file[,c(2,rem_cols)])
for (j in 2:(ncol(ratio_file_mCherry)-1)){
  col1 <- as.matrix(ratio_file_mCherry[,c(j)])
  col_back <- as.matrix(ratio_file_mCherry[,c(ncol(ratio_file_mCherry))])
  ratio_file_mCherry[,c(j)]<-col1-col_back
}
ratio_file_mCherry[,c(ncol(ratio_file_mCherry))]<-NULL

# The following block corrects bleach of mCherry fluorescence (empirically found to occur during the recordings)

for (i in 2:ncol(ratio_file_mCherry)){
  full_data <- as.data.frame(ratio_file_mCherry[c(1:which(ratio_file_GFP$X1==1)),c(1,i)])
  full_data2 <- as.data.frame(ratio_file_mCherry[,c(1,i)])
  colnames(full_data)<-c("X","Y")
  colnames(full_data2)<-c("X","Y")
  nlsfit <- nls(formula = Y ~ Y0*exp(k * X), data = full_data, start = list(k = 0.01,Y0=max(full_data$Y)))
  nlsfit_vals <- as.data.frame(tidy(nlsfit))
  Y0 <- as.numeric(nlsfit_vals[2,2])
  k <- as.numeric(nlsfit_vals[1,2])
  full_data2 <- full_data2 %>% mutate (Y2=full_data2$Y+(Y0-Y0*exp(k*full_data2$X)))
  ratio_file_mCherry[,c(i)]<-full_data2$Y2
}
ratio_file_mCherry[,c(1)]<-NULL

# The following block calculates the GCaMP/tdTomato ratio

final_ratios <- as.data.frame(matrix(0,nrow=nrow(ratio_file_GFP),ncol=ncol(ratio_file_GFP)))
final_ratios[,c(1,2)]<-ratio_file_GFP[,c(1,2)]
colnames(final_ratios)<-colnames(ratio_file_GFP)
for (k in 3:ncol(final_ratios)){
  GFP<-ratio_file_GFP[,c(k)]
  mCherry<-ratio_file_mCherry[,c(k-2)] 
  final_ratios[,c(k)]<-GFP/mCherry
}

# The following block sets the baseline (from beginning of experiment to event 1) to 1.

for (i in 3:ncol(final_ratios)){
  col1 <- as.data.frame(final_ratios[,c(i)])
  col2 <- col1[c(1:(which(final_ratios$X1==1))),]
  norm_val <- mean(col2)
  final_ratios[,c(i)]<-final_ratios[,c(i)]/norm_val
}

event_number <- 1 # Event number that marks the starting time for AUC calculation.
duration_expt <- 120 # Duration (in seconds) for which AUC is calculated. This would be after the event specified in the previous line.

area_list <- as.data.frame(matrix(0,nrow=0,ncol=1))
for (k in 3:ncol(final_ratios)){
  col1 <- as.data.frame(final_ratios[,c(1,2,k)])
  colnames(col1)<-c("event","time","ratio")
  col1 <- col1 %>% mutate (base=1)
  col1[c(which(col1$ratio<1)),3]<-1 # This line is meant to remove negative AUC values.
  time1 <- col1[c(which(col1$event==event_number)),c(2)]
  time2 <- min(col1[c(which(col1$time>(time1+duration_expt))),c(2)])
  area1 <- auc(col1$time,col1$ratio,from = time1,to=time2)
  area2 <- auc(col1$time,col1$base,from = time1,to=time2)
  area <- as.data.frame(area1-area2)
  colnames(area)<-colnames(area_list)
  area_list <- as.data.frame(rbind(area_list,area))
}

write.csv(area_list,"~/AUC_results.csv") # The name of the output file can be modified as needed
