#This script will create statistics of geneclass assignment 
#to determine how well this analysis performed

######################
#     Libraries      #
######################

###############################
#     Load in data files      #
###############################

setwd("../../QUAC_data_files")

###########################
#     Load Data Files     #
###########################

#load in file with sample pop info 
QUAC_sample_pops <- read.csv("QUAC_geneclass/QUAC_sample_pop.csv")

#load in relate red df 
QUAC_relate_red_df <- read.csv("QUAC_data_frames/Relate_Red/QUAC_garden_relate_red_df.csv")

#limit sample pops df by relate_red 
QUAC_sample_pops_relate_red <- QUAC_sample_pops[QUAC_sample_pops$Ind %in% QUAC_relate_red_df$Ind,]

###################################
#     Assess assignment test      #
###################################
setwd("QUAC_geneclass")
#load in all output files 
QUAC_assign_output <- list.files(pattern = "_output.csv")

#create store list 
QUAC_assign_list <- list()

#create final matrix to store results
QUAC_assign_success <- matrix(nrow = length(QUAC_assign_output), ncol = 1)

#loop to load in assignment output files 
for(assign in 1:length(QUAC_assign_output)){
  
  #load in assignment files 
  QUAC_assign_list[[assign]] <- read.csv(paste0(QUAC_assign_output[[assign]]))
  
  #add sample population depending on if they are reduced or not 
  if(assign == 1|assign == 2){
    
    QUAC_assign_list[[assign]] <- cbind(QUAC_assign_list[[assign]], QUAC_sample_pops$Sample_Pop) 
    
  }else{
    
    QUAC_assign_list[[assign]] <- cbind(QUAC_assign_list[[assign]], QUAC_sample_pops_relate_red$Sample_Pop) 
  }
  
  #now determine correct assignment - if the assignment testing correctly assigned individuals
  #add a "Y" if no, add a "N", if no sample pop, say "none"
  
  for(ind in 1:length(QUAC_assign_list[[assign]][,1])){
  assign_test <- QUAC_assign_list[[assign]][ind,13] == QUAC_assign_list[[assign]][ind,3]
  
  if(assign_test == TRUE){
    
    QUAC_assign_list[[assign]][ind,14] <- "Y"
    
  }else{
    
    QUAC_assign_list[[assign]][ind,14] <- "N"
    
    }
    
  }
  
  #limit data frame by individuals with no source pop info 
  QUAC_lm_df_1 <- QUAC_assign_list[[assign]][QUAC_assign_list[[assign]][,13] == "NONE",]
  
  #remove these from the data frame 
  QUAC_lm_df2 <- QUAC_assign_list[[assign]][!QUAC_assign_list[[assign]][,13] %in% QUAC_lm_df_1[,13],]
  
  #now calculate the final result 
  QUAC_assign_success[assign,1] <- signif(length(QUAC_lm_df2[QUAC_lm_df2[,14] == "Y",][,1])/length(QUAC_lm_df2[,1])*100, 3)
 
  #name rows 
  rownames(QUAC_assign_success) <- gsub("\\..*", "",QUAC_assign_output)
}

#write out final data file 
setwd("../../QUAC_analyses/Results/Clustering")

write.csv(QUAC_assign_success, "QUAC_assign_success.csv")



