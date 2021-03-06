# ANALYSIS OF PERCEVALHR DATA

library(tidyverse)
library(readr)
library(MESS)

# Data file to be used should be a CSV formatted as indicated here. 
# Column 1 should contain the numbers 1, 2, 3.. that mark the events such as oligomycin A addition. This column should NOT have a label.
# Column 2 should contain time in seconds. This column can have any label.
# Column 3 onward should contain the PercevalHR ratios. These columns can have any label.
# Values from the beginning to the first addition event (1 in column 1) are considered the baseline and will be set to 1

ratio_file <- read_csv("~/ratio.csv") # The name of the input file can be modified as needed
ratio_file2 <- ratio_file

# The following block sets the baseline to 1.

for (i in 3:ncol(ratio_file2)){
  col1 <- as.data.frame(ratio_file2[,c(i)])
  col2 <- col1[c(1:(which(ratio_file2$X1==1)-1)),]
  norm_val <- mean(col2)
  ratio_file2[,c(i)]<-ratio_file2[,c(i)]/norm_val
}
ratio_file3 <- as.data.frame(ratio_file2[-c(1:(which(ratio_file2$X1==1)-1)),])
ratio_file3[,c(2)]<-ratio_file3[,c(2)]-ratio_file3[c(1),c(2)]

area_list <- as.data.frame(matrix(0,nrow=0,ncol=1))

duration_of_expt <- 120 # this is an important parameter that establishes the duration, after the addition event, over which AUC will be calciulated 

# The following block determines AUC for each column with PercevalHR ratios.

for (k in 3:ncol(ratio_file3)){
  col1 <- as.data.frame(ratio_file3[,c(1,2,k)])
  colnames(col1)<-c("adds","time","ratio")
  col1 <- col1 %>% mutate (base=1)
  time1 <- as.numeric(col1[c(which(col1$adds==1)),2])
  time2 <- as.numeric(col1[c(min(which(col1$time>(time1+duration_of_expt)))-1),2])
  area1 <- auc(col1$time,col1$ratio,from = time1,to=time2)
  area2 <- auc(col1$time,col1$base,from = time1,to=time2)
  area <- as.data.frame(area1-area2)
  colnames(area)<-colnames(area_list)
  area_list <- as.data.frame(rbind(area_list,area))
}

write.csv(area_list,"~/AUC_results.csv") # The name of the output file can be modified as needed
