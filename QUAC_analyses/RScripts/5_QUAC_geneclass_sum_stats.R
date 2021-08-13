#This script will compare the statistics of geneclass assignment 

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


QUAC_assignment_output <- list.files(path = "G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\GeneClass", 
                                     pattern = "output.csv$")
##list of data frames 
QUAC_assignment_output_df <- list()

##data frame to store all of the results 
QUAC_classification_percent <- matrix(nrow = 3, ncol = 4)

##list of data frames 
QUAC_classification_samp_pop <- matrix(nrow = length(rownames(QUAC_assignment_output_df[[1]])), ncol = 3)

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
  QUAC_classification_str_samp_pop[i,] <- length(QUAC_assignment_df[(QUAC_classification_samp_pop[,i] == "TRUE"),][,1])/length(QUAC_classification_samp_pop[,1])*100
}

