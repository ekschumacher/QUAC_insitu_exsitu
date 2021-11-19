#This script will create statistics of geneclass assignment 
#to determine how well this analysis performed

##########################
######## Libraries #######
##########################

library(adegenet)

#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_insitu_exsitu\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_insitu_exsitu\\QUAC_analyses\\Results"

###############################################
####### Read in classification data frame #####
###############################################
setwd(QUAC_data_files)

##read in data files 
QUAC_assignment_output_df <- read.csv(paste0(QUAC_analysis_results, "\\Clustering\\Assignment_red.csv"))

##now figure out # of times we have a yes 
QUAC_assignment_correct <- na.omit(QUAC_assignment_output_df[QUAC_assignment_output_df$Correct_Assignment == "Y",])

##now figure out % 
length(QUAC_assignment_correct[,1])/length(QUAC_assignment_output_df[,1])*100
##assigned to the correct population 25.3% of the time 

##now try with structure pops 
QUAC_assignment_correct_str_df <- read.csv(paste0(QUAC_analysis_results, "\\Clustering\\assignment_red_3pop.csv"))

##now figure out number of yes 
QUAC_assignment_correct_str <- na.omit(QUAC_assignment_correct_str_df[QUAC_assignment_correct_str_df$Correct_Assignment == "Y",])

##now calculate the percent 
length(QUAC_assignment_correct_str[,1])/length(QUAC_assignment_correct_str_df[,1])*100

##load in output data frames from Geneclass
setwd(paste0(QUAC_analysis_results, "\\Clustering"))
QUAC_assignment_output <- list.files(pattern = "assignment_")
##list of data frames 
QUAC_assignment_output_df <- list()

##data frame to store all of the results 
QUAC_classification_percent <- matrix(nrow = 2, ncol = 4)

##loops to calculate % classification
for(i in 1:length(QUAC_assignment_output)){
  ##data frame loading 
  QUAC_assignment_output_df[[i]] <- read.csv(QUAC_assignment_output[[i]])
  
  ##fill in data frame for percent classification 
  QUAC_classification_percent[i,1] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,6] >= 25,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
  
  QUAC_classification_percent[i,2] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,6] >= 50,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
  
  ##75% of the time 
  QUAC_classification_percent[i,3] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,6] >= 75,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
  
  ##95% of the time 
  QUAC_classification_percent[i,4] <- length(rownames(QUAC_assignment_output_df[[i]][QUAC_assignment_output_df[[i]][,6] >= 95,]))/(length(rownames(QUAC_assignment_output_df[[i]])))*100
 
}

#name rows and columns 
rownames(QUAC_classification_percent) <- c("All Wild Pops", "Three Str Pops")
colnames(QUAC_classification_percent) <- c("Assign_25%","Assign_50%","Assign_75%","Assign_95%")

##write out data frames
write.csv(QUAC_classification_percent, paste0(QUAC_analysis_results, "\\Clustering\\percent_assignment.csv"))
