##############This script details the preliminary analyses
#############conducted on Quercus acerifolia 
############(referred to as QUAC for brevity here)
###########For the project looking at genetic diversity capture 
##########in garden populations compared to wild populations

##########################
######## Libraries #######
##########################

library(diveRsity)
library(adegenet)
#library(tidyr)
#library(poppr)
#library(hierfstat)

#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

##############################################
############# Convert to Genind ##############
##############################################

setwd(QUAC_data_files)

##create a list of arlequin files
QUAC_arp_list <- list.files(pattern = ".arp$")

##create a list to convert to genind files
for(a in 1:length(QUAC_arp_list)){
  
  arp2gen(QUAC_arp_list[[a]])
  
}

##now read in list of genind files 
QUAC_gen_list <- list.files(pattern = ".gen$")

##finally, create a list of genind files here 
QUAC_genind_list <- list()

##loop over to load in genind files 
for(g in 1:length(QUAC_gen_list)){
  
  QUAC_genind_list[[g]] <- read.genepop(QUAC_gen_list[[a]], ncode = 3)
  
}







