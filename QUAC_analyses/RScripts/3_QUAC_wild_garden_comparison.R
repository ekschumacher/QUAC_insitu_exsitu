##########################
######## Libraries #######
##########################

library(adegenet)
library(diveRsity)

#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

################################
########## Load files ##########
################################
setwd(QUAC_data_files)

##convert to a genind 
arp2gen(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_garden_wild_clean.arp"))

##load in genind for QUAC
QUAC_wild_garden_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_garden_wild_clean.gen"), ncode = 3)

##load in data frame 
QUAC_wild_garden_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\QUAC_wild_garden_df.csv"))

##rename individuals in the genind file
rownames(QUAC_wild_garden_gen@tab) <- QUAC_wild_garden_df$Ind

##rename population names 
QUAC_popnames <- unique(QUAC_wild_garden_df$Pop)

levels(QUAC_wild_garden_gen@pop) <- QUAC_popnames



