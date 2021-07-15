#These are the next steps in calculating genetic diversity 
#First, we tested linkage disequilibrium, null alleles, 
#and Hardy Weinberg equilibrium.
#Next, a data frame including expected heterozygosity, allelic richness,
#number of alleles, mean longtiude and latitude by population, and 
#individual numbers. This table is included in full in the supplemental text 
#of this manuscript.

##########################
######## Libraries #######
##########################

library(diveRsity)
#library(adegenet)
#library(poppr)
#library(Demerelate)

#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

#######################################
########## Load genind files ##########
#######################################

setwd(QUAC_data_files)

arp2gen(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_allpop_clean.arp"))


