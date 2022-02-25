##This script details an assessment of the number of individuals
#sourced from the same maternal plant in botanic garden 
#populations.

#########################
#        Libraries      #
#########################

library(tidyverse)
library(Demerelate)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("../../QUAC_data_files")

#load in accession records 
QUAC_garden_accessions <- read.csv("QUAC_data_frames/QUAC_garden_accessions.csv")

###################################
#     Maternal plant records      #
###################################
#first reduce accession records to just leave the "mother plant" indication
QUAC_garden_accessions <- QUAC_garden_accessions %>% 
                            mutate(accession_simple = gsub("\\*.*","",QUAC_garden_accessions$Accession))

#code to calculate unique maternal lines 
QUAC_total_mom <- length(unique(QUAC_garden_accessions$accession_simple))

#create botanic garden list 
QUAC_bg_list <- unique(QUAC_garden_accessions$Pop)

#create a data frame to store maternal lines by BG
QUAC_bg_maternal <- matrix(nrow = length(QUAC_bg_list), ncol = 3)

#load in QUAC relatedness analysis ready DF 
QUAC_df <- read.csv("QUAC_data_frames/Garden_Wild/QUAC_clean_df.csv")
#reduce to just garden 
QUAC_rel_ready_df <- QUAC_df[QUAC_df$Garden_Wild == "Garden",]
#also remove weird column 
QUAC_rel_ready_df <- QUAC_rel_ready_df[,-1]
QUAC_rel <- Demerelate(QUAC_rel_ready_df[,-2], object = T, value = "loiselle")
##reduce individuals based on relatedness
#extract just the names of highly related inds (half-sib or more)
QUAC_halfsib_names <- names(which(unlist(QUAC_rel$Empirical_Relatedness) > 0.25))

##remove extra characters
QUAC_halfsib_names_cleanfront <-  gsub("^.*\\.","", QUAC_halfsib_names)

##remove extra name
QUAC_halfsib_names_cleanback <- gsub("^.*\\_","", QUAC_halfsib_names_cleanfront)

##create list of individuals to remove 
QUAC_ind_red_list <- unique(QUAC_halfsib_names_cleanback)

##reduce data frame for related individuals 
QUAC_relate_red_df <- QUAC_rel_ready_df[!QUAC_rel_ready_df[,1] %in% QUAC_ind_red_list,]


#code calculate maternal lines by garden 
for(garden in 1:length(QUAC_bg_list)){
  

  #count number of individuals per botanic garden
  QUAC_bg_maternal[garden,1] <- length(QUAC_garden_accessions[QUAC_garden_accessions$Pop == paste0(QUAC_bg_list[[garden]]),][,1])
  
  #count unique maternal lines for each botanic garden
  QUAC_bg_maternal[garden,2] <- length(unique(QUAC_garden_accessions[QUAC_garden_accessions$Pop == paste0(QUAC_bg_list[[garden]]),]$accession_simple))
  
  #count the number of individuals with genetic similarities 
  QUAC_bg_maternal[garden,3] <- length(QUAC_relate_red_df[QUAC_relate_red_df$Pop == QUAC_bg_list[[garden]],][,1])
  
}


#name rows and columns
rownames(QUAC_bg_maternal) <- QUAC_bg_list
colnames(QUAC_bg_maternal) <- c("Total_Ind", "Unique_Maternal_Lines", "Genetically_Distinct_Ind")

#write out 
write.csv(QUAC_bg_maternal, "../QUAC_analyses/Results/Sum_Stats/QUAC_bg_maternal.csv")

