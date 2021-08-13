#This script will create statistics of geneclass assignment 
#to determine how well this analysis performed

##########################
######## Libraries #######
##########################

library(adegenet)

#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

###############################################
####### Read in classification data frame #####
###############################################
setwd(QUAC_data_files)
##load in output data frames from Geneclass
QUAC_assignment_output <- list.files(path = "QUAC_geneclass", pattern = "output.csv$")
##list of data frames 
QUAC_assignment_output_df <- list()

##data frame to store all of the results 
QUAC_classification_percent <- matrix(nrow = 3, ncol = 4)

##list of data frames 
QUAC_classification_samp_pop <- matrix(nrow = 146, ncol = 3)

##data frame to store % correct assignment for structured assignment 
QUAC_classification_str_samp_pop <- matrix(nrow = 3, ncol = 1)

##loops to calculate % classification
for(i in 1:length(QUAC_assignment_output)){
  ##data frame loading 
  QUAC_assignment_output_df[[i]] <- read.csv(paste0("G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\GeneClass\\", QUAC_assignment_output[[i]]))
  
  ##fill in data frame for percent classification 
  QUAC_classification_percent[i,1] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,5] >= 25,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
  
  QUAC_classification_percent[i,2] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,5] >= 50,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
  
  ##75% of the time 
  QUAC_classification_percent[i,3] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,5] >= 75,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
  
  ##95% of the time 
  QUAC_classification_percent[i,4] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,5] >= 95,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
  
  ##now create a data frame to compare if part of a word is contained in the data frame 
  for(j in 1:146) QUAC_classification_samp_pop[j,i] <- grepl(QUAC_assignment_output_df[[i]][j,3], QUAC_assignment_output_df[[i]][j,4], fixed = TRUE)
  
  ##calculate the % correct assignment rate in the first column 
  QUAC_classification_str_samp_pop[i,] <- length(QUAC_classification_samp_pop[(QUAC_classification_samp_pop[,i] == "TRUE"),][,i])/length(QUAC_classification_samp_pop[,i])*100

}

#name rows and columns 
rownames(QUAC_classification_percent) <- c("WildPops","SixPops","ThreePops")
colnames(QUAC_classification_percent) <- c("Assign_25%","Assign_50%","Assign_75%","Assign_95%")
rownames(QUAC_classification_str_samp_pop) <- c("WildPops","SixPops","ThreePops")
colnames(QUAC_classification_str_samp_pop) <- "Percent_Classified_Correctly"

##write out data frames
write.csv(QUAC_classification_percent, paste0(QUAC_analysis_results, "\\Clustering\\Assignment_percentage.csv"))
write.csv(QUAC_classification_str_samp_pop, paste0(QUAC_analysis_results, "\\Clustering\\pop_correct_assign.csv"))

