##This script details the analyses run on Quercus acerifolia (referred to as 
#QUAC for brevity here) genotype files to prepare them for genetic diversity and 
#structure analyses. When files are referred to as "clean" that means individuals 
#that are clones and individuals with too much missing data have been removed. 
#When files and objects are titled "red" that means they have been reduced
#for relatedness (25% or more related individuals are reduced to one individual
#per phenotype)

#########################
#        Libraries      #
#########################

library(diveRsity)
library(adegenet)
library(poppr)
library(Demerelate)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("../../QUAC_data_files")

##convert QUAC all pop arlequin file into genepop file if needed
#arp2gen("QUAC_adegenet_files/QUAC_allpop.arp")

##now read in genepop file as a genind for adegenet 
QUAC_gen <- read.genepop("QUAC_adegenet_files/QUAC_allpop.gen", ncode = 3)

##load relatedness data frame for relatedness analysis 
QUAC_df <- read.csv("QUAC_data_frames/QUAC_allpop.csv")
QUAC_df <- QUAC_df[,-1]

##rename individuals in the genind object 
rownames(QUAC_gen@tab) <- QUAC_df[,1]

#load in relatedness function 
source("../QUAC_analyses/RScripts/relatedness_analyses.R")

############################################################
#    Remove Clones and Individuals with Missing Data       #
############################################################
##clone check 
#convert to genelcone object
QUAC_geneclone <- as.genclone(QUAC_gen)
#identify multi-locus genotypes (non-clones)
QUAC_gen_id <- mlg.id(QUAC_geneclone)
#function to pull out all clones into a list
QUAC_clone_index <- which(sapply(QUAC_gen_id, function(x) length(x)>1))
#create a list of ID with clones
clones_ID <- list()
for(clones in 1:length(QUAC_clone_index)) clones_ID[[clones]] <- QUAC_gen_id[[QUAC_clone_index[[clones]]]]
#now remove clones from the matrix 
QUAC_noclones <- clonecorrect(QUAC_geneclone)
#convert back to a genind object
QUAC_genind_nocl <- genclone2genind(QUAC_noclones) ##left with 456 individuals, there were 7 clones 

##remove individuals with too much missing data 
QUAC_genind_nomd <- missingno(QUAC_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) ##there were 7 individuals with too much missing data, left with 449 individuals

#output to genalex csv if needed 
genind2genalex(QUAC_genind_nomd, file="QUAC_data_frames/Garden_Wild/QUAC_garden_wild_clean_genalex.csv", 
              overwrite = TRUE)

##remove the individuals from the relatedness df that have been removed for missing data and clones
QUAC_clean_df <- QUAC_df[QUAC_df[,1] %in% rownames(QUAC_genind_nomd@tab),]

##write out data frame 
write.csv(QUAC_clean_df, "QUAC_data_frames/Garden_Wild/QUAC_clean_df.csv")

#################################
#      Relatedness Analysis     #
#################################
##create pop type list 
pop_type_list <- c("Garden", "Wild")

##loop to calculate relatedness for garden and wild individuals and reduce data files for relatedness
for(pop_type in 1:length(pop_type_list)){
  
  #limit data frame by pop_type
  QUAC_df_temp <- QUAC_clean_df[QUAC_clean_df$Garden_Wild == paste0(pop_type_list[[pop_type]]),]
  
  #limit genind object by pop_type 
  QUAC_gen_temp <- QUAC_genind_nomd[rownames(QUAC_genind_nomd@tab) %in% QUAC_df_temp[,1],]
  
  #calculate relatedness between all individuals for either wild or garden df
  QUAC_rel_df <- halfsib_relatedness_reduction_loiselle(QUAC_df_temp[,-2])
  
  #write csv of the reduced data frame 
  write.csv(QUAC_rel_df, paste0("QUAC_data_frames/Relate_Red/QUAC_", pop_type_list[[pop_type]], "_relate_red_df.csv"))
  
  #now limit genind object by relatedness 
  QUAC_relate_red_gen <- QUAC_gen_temp[rownames(QUAC_gen_temp@tab) %in% QUAC_rel_df$Ind,]
  
  #write out genalex file
  genind2genalex(QUAC_relate_red_gen, paste0("QUAC_data_frames/Relate_Red/QUAC_", pop_type_list[[pop_type]], 
                                    "_relate_red_genalex.csv"), overwrite = TRUE)
  
}

