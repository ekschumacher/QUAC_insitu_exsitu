#This script will create statistics of geneclass assignment 
#to determine how well this analysis performed

##########################
######## Libraries #######
##########################

library(adegenet)

###############################################
####### Read in classification data frame #####
###############################################
setwd("../../QUAC_data_files")

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

#########################################
######## Assessing performance ##########
#########################################
setwd("../QUAC_analyses/Results/Clustering")

QUAC_allgarden_output <- read.csv("G:/Shared drives/Emily_Schumacher/QUAC_analyses/QUAC_geneclass_12_16_21/QUAC_garden_wild_assignment_results.csv")
QUAC_red_relate_garden_output <- read.csv("QUAC_red_relate_garden_output.csv")
QUAC_allgarden_red_relate_wild <- read.csv("QUAC_wild_red_relate_allgarden_output_geneclass.csv")

##create a new column
QUAC_allgarden_output[,14] <- NA

for(n in 1:length(QUAC_allgarden_output[,1])){
  
  h <- QUAC_allgarden_output[n,3] == QUAC_allgarden_output[n,4]
  
  if(h == TRUE){
    
    QUAC_allgarden_output[n,14] <- "Y"
    
  }else{
    
    QUAC_allgarden_output[n,14] <- "N"
    
  }
  
}

##then limit by the individuals with source populations
QUAC_allgarden_df2 <- QUAC_allgarden_output[QUAC_allgarden_output$Sample_pop == "NONE",]

##then limit the df by those rows
QUAC_allgarden_output <- QUAC_allgarden_output[!QUAC_allgarden_output$Ind %in% QUAC_allgarden_df2$Ind,]

##
QUAC_ass_test_success_rate <- length(QUAC_allgarden_output[QUAC_allgarden_output[,14] == "Y",][,1])/length(QUAC_allgarden_output[,1])*100
##about 25%

##create a new column
QUAC_red_relate_garden_output[,17] <- NA

##run the loop for testing with individuals reduced 
for(n in 1:length(QUAC_red_relate_garden_output[,1])){
  
  h <- QUAC_red_relate_garden_output[n,3] == QUAC_red_relate_garden_output[n,4]
  
  if(h == TRUE){
    
    QUAC_red_relate_garden_output[n,17] <- "Y"
    
  }else{
    
    QUAC_red_relate_garden_output[n,17] <- "N"
    
  }
  
}
##limit by the individuals with no source info 
QUAC_red_relate_df2 <- QUAC_red_relate_garden_output[QUAC_red_relate_garden_output$Sample_Pop == "NONE",]
##then remove those from the output df
QUAC_red_relate_garden_output <- QUAC_red_relate_garden_output[!QUAC_red_relate_garden_output$Ind %in% QUAC_red_relate_df2$Ind,]

QUAC_ass_test_success_rate_red_relate_garden <- length(QUAC_red_relate_garden_output[QUAC_red_relate_garden_output[,17] == "Y",][,1])/length(QUAC_red_relate_garden_output[,1])*100

##try assignment test with wild individuals reduced 
QUAC_allgarden_red_relate_wild[,17] <- NA

##run the loop for testing with individuals reduced 
for(n in 1:length(QUAC_allgarden_red_relate_wild[,1])){
  
  h <- QUAC_allgarden_red_relate_wild[n,3] == QUAC_allgarden_red_relate_wild[n,4]
  
  if(h == TRUE){
    
    QUAC_allgarden_red_relate_wild[n,17] <- "Y"
    
  }else{
    
    QUAC_allgarden_red_relate_wild[n,17] <- "N"
    
  }
  
}

QUAC_ass_test_success_rate_red_relate_wild_allgarden <- length(QUAC_allgarden_red_relate_wild[QUAC_allgarden_red_relate_wild[,17] == "Y",][,1])/length(QUAC_allgarden_red_relate_wild[,1])*100

##load in data file
QUAC_red_relate_wild_woKessler_allgarden <- read.csv("QUAC_wild_red_relate_woKessler_allgarden_output_geneclass.csv")

##now create a new column for the classification test
QUAC_red_relate_wild_woKessler_allgarden[,17] <- NA

##run the loop for testing with individuals reduced 
for(n in 1:length(QUAC_red_relate_wild_woKessler_allgarden[,1])){
  
  h <- QUAC_red_relate_wild_woKessler_allgarden[n,3] == QUAC_red_relate_wild_woKessler_allgarden[n,4]
  
  if(h == TRUE){
    
    QUAC_red_relate_wild_woKessler_allgarden[n,17] <- "Y"
    
  }else{
    
    QUAC_red_relate_wild_woKessler_allgarden[n,17] <- "N"
    
  }
  
}

QUAC_df2 <- QUAC_red_relate_wild_woKessler_allgarden[QUAC_red_relate_wild_woKessler_allgarden$Sample_Pop == "NONE",]
##reduc
QUAC_red_relate_wild_woKessler_allgarden <- QUAC_red_relate_wild_woKessler_allgarden[!QUAC_red_relate_wild_woKessler_allgarden$Ind %in% QUAC_df2$Ind,]

QUAC_ass_test_success_rate_red_wild_woKESSLER <- length(QUAC_red_relate_wild_woKessler_allgarden[QUAC_red_relate_wild_woKessler_allgarden[,17] == "Y",][,1])/length(QUAC_red_relate_wild_woKessler_allgarden[,1])*100
##38.1%

##now create a new column for the classification test
QUAC_allgarden_red_relate_wild[,17] <- NA

##run the loop for testing with individuals reduced 
for(n in 1:length(QUAC_allgarden_red_relate_wild[,1])){
  
  h <- QUAC_allgarden_red_relate_wild[n,3] == QUAC_allgarden_red_relate_wild[n,4]
  
  if(h == TRUE){
    
    QUAC_allgarden_red_relate_wild[n,17] <- "Y"
    
  }else{
    
    QUAC_allgarden_red_relate_wild[n,17] <- "N"
    
  }
  
}
##remove individuals without source info
QUAC_allgard_df2 <- QUAC_allgarden_red_relate_wild[QUAC_allgarden_red_relate_wild$Sample_Pop == "NONE",]
##reduc
QUAC_allgarden_red_relate_wild <- QUAC_allgarden_red_relate_wild[!QUAC_allgarden_red_relate_wild$Ind %in% QUAC_allgard_df2$Ind,]

QUAC_ass_test_success_allgard <- length(QUAC_allgarden_red_relate_wild[QUAC_allgarden_red_relate_wild[,17] == "Y",][,1])/length(QUAC_allgarden_red_relate_wild[,1])*100

##read in the all wild, reduced garden individuals
QUAC_wild_red_garden <- read.csv("QUAC_allwild_red_garden_geneclass.csv")

##now create a new column for the classification test
QUAC_wild_red_garden[,13] <- NA

##run the loop for testing with individuals reduced 
for(n in 1:length(QUAC_wild_red_garden[,1])){
  
  h <- QUAC_wild_red_garden[n,2] == QUAC_wild_red_garden[n,3]
  
  if(h == TRUE){
    
    QUAC_wild_red_garden[n,13] <- "Y"
    
  }else{
    
    QUAC_wild_red_garden[n,13] <- "N"
    
  }
  
}
##now remove anything with no sample pop info 
QUAC_wild_red_gard_df2 <- QUAC_wild_red_garden[QUAC_wild_red_garden$Sample_Pop == "NONE",]
##limit data frame by this 
QUAC_wild_red_garden <- QUAC_wild_red_garden[!QUAC_wild_red_garden$Ind %in% QUAC_wild_red_gard_df2$Ind,]

##now calculate the % correct loading
QUAC_ass_test_success_redgarden_allwild <- length(QUAC_wild_red_garden[QUAC_wild_red_garden[,13] == "Y",][,1])/length(QUAC_wild_red_garden[,1])*100
