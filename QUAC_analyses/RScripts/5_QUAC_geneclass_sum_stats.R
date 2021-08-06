#This script will compare the statistics of geneclass assignment 

##########################
######## Libraries #######
##########################



#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

###############################################
####### Read in classification data frame #####
###############################################

QUAC_assignment_output <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\QUAC_assignment_output.csv"))

##calculate the percent of the individuals assignment 25% of the time 
length(rownames(QUAC_assignment_output[QUAC_assignment_output$Percent_1 >= 25,]))/(length(rownames(QUAC_assignment_output)))*100

##50% of the time 
length(rownames(QUAC_assignment_output[QUAC_assignment_output$Percent_1 >= 50,]))/(length(rownames(QUAC_assignment_output)))*100

##75% of the time 
length(rownames(QUAC_assignment_output[QUAC_assignment_output$Percent_1 >= 75,]))/(length(rownames(QUAC_assignment_output)))*100

##95% of the time 
length(rownames(QUAC_assignment_output[QUAC_assignment_output$Percent_1 >= 95,]))/(length(rownames(QUAC_assignment_output)))*100

###
QUAC_assignment_df <- matrix(nrow = length(QUAC_assignment_output$Ind), ncol = 2)

##
for(i in 1:length(QUAC_assignment_output$Ind)){
  
  QUAC_assignment_df[i,1] <- identical(QUAC_assignment_output[i,3],QUAC_assignment_output[i,4])
  
  QUAC_assignment_df[i,2] <- identical(QUAC_assignment_output[i,3],QUAC_assignment_output[i,6])
}

length(QUAC_assignment_df[(QUAC_assignment_df[,1] == "TRUE"),][,1])/length(QUAC_assignment_df[,1])*100

length(QUAC_assignment_df[(QUAC_assignment_df[,1] == "TRUE")|(QUAC_assignment_df[,2] == "TRUE"),][,1])/length(QUAC_assignment_df[,1])*100

QUAC_S1_df <- QUAC_assignment_output[QUAC_assignment_output$Sampled.From == c("Kessler","Magazine","SL","Porter"),]
QUAC_S1_df$Assign <- "S1"

QUAC_S2_df <- QUAC_assignment_output[QUAC_assignment_output$Sampled.From == "Pryor",]
QUAC_S2_df$Assign <- "S2"
