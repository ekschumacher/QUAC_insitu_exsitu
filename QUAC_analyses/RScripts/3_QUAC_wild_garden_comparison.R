##This script details the analysis of actually comparing genetic 
#diversity levels between garden and wild Q. acerifolia populations. 
#We first calculated diversity levels throughout all garden  
#and wild populations (indicated by allelic richness and 
#expected heterozygosity) and then ran a t-test to assess significance. 
#In this script we use "cleaned" genetic files because they have been
#cleaned for clones and individuals with too much misisng data 
#(25% or more)

#####################
#     Libraries     #
#####################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)

#######################
#     Load files      #
#######################
setwd("../../QUAC_data_files")

##convert to a genepop file, if necessary
#arp2gen(paste0(QUAC_data_files, "/QUAC_genind/QUAC_garden_wild_clean.arp"))

#load in genepop file as a genind object 
QUAC_garden_wild_gen <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_garden_wild_clean.gen", ncode = 3)

#load in data frame 
QUAC_garden_wild_df <- read.csv("QUAC_data_frames/Garden_Wild/QUAC_clean_df.csv")[,-1]

#rename individuals in genind object
rownames(QUAC_garden_wild_gen@tab) <- QUAC_garden_wild_df$Ind

##rename population names in the genind object
#create pop name list
QUAC_popnames <- unique(QUAC_garden_wild_df$Pop)
#rename populations in the genind object
levels(QUAC_garden_wild_gen@pop) <- QUAC_popnames

#load in function to calculate allele frequency categories
source("../QUAC_analyses/RScripts/Fa_sample_funcs.R")

#create functions to run code 
colMax <- function(data) sapply(data, max, na.rm = TRUE)

#################################################
#     Comparing wild and garden populations     #
#################################################
#set working directory 
setwd("../QUAC_analyses/Results/Wild_Garden_Comparison")

#allrich list 
QUAC_allrich <- list()
QUAC_hexp <- list()

#genetic diversity data frame 
QUAC_gendiv_df <- matrix(nrow = 3, ncol = 2)

##name rows and columns 
rownames(QUAC_gendiv_df) <- c("Mean_Garden", "Mean_Wild", "p-value")
colnames(QUAC_gendiv_df) <- c("All_Rich","Hexp")

#loop over both garden and wild individuals to calculate the differences in allelic richness and hexp
for(pop_type in 1:length(pop_type_list)){
  
  #create data frame of the pop type - garden or wild 
  QUAC_df <- QUAC_garden_wild_df[QUAC_garden_wild_df$Garden_Wild == paste0(pop_type_list[[pop_type]]),]
  
  #then subset genind by pop type - garden of wild 
  QUAC_gen <- QUAC_garden_wild_gen[rownames(QUAC_garden_wild_gen@tab) %in% QUAC_df$Ind,]
  
  ####Allelic richness calculations 
  QUAC_allrich[[pop_type]] <- allelic.richness(QUAC_gen)$Ar
  
  if(pop_type == 1){
    QUAC_allrich[[1]]$pop_type <- "Garden"
  }else{
    QUAC_allrich[[2]]$pop_type <- "Wild"
  }
  #combine to form an allelic richness data frame 
  QUAC_allrich_df <- rbind(QUAC_allrich[[1]][,c(1,3)], QUAC_allrich[[2]][,c(1,3)])
  #add allelic richness results to the data frame 
  QUAC_gendiv_df[pop_type,1] <- mean(QUAC_allrich_df[QUAC_allrich_df[,2] == paste0(pop_type_list[[pop_type]]),][,1])
  QUAC_gendiv_df[3,1] <- kruskal.test(QUAC_allrich_df[,1]~QUAC_allrich_df[,2])[3]$p.value

  ###Expected heterozygosity 
  QUAC_hexp[[pop_type]] <- data.frame(summary(QUAC_gen)[7]$Hexp)
  
  #name the pops 
  if(pop_type == 1){
    QUAC_hexp[[1]]$pop_type <- "Garden"
  }else{
    QUAC_hexp[[2]]$pop_type <- "Wild"
  }
  #combine all of the heterozygosity calculations 
  QUAC_hexp_df <- rbind(QUAC_hexp[[1]][,c(1:2)], QUAC_hexp[[2]][,c(1:2)])
  #add allelic richness results to the data frame 
  QUAC_gendiv_df[pop_type,2] <- mean(QUAC_hexp_df[QUAC_hexp_df[,2] == paste0(pop_type_list[[pop_type]]),][,1])
  QUAC_gendiv_df[3,2] <- kruskal.test(QUAC_hexp_df[,1]~QUAC_hexp_df[,2])[3]$p.value
  
  #write out a csv 
  write.csv(QUAC_gendiv_df, "QUAC_gendiv_df.csv")
}

#################################
#     Allelic capture code      #
#################################
##list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

##seppop genind - we only want to use the wild pops to calculate the alleles existing for capture
QUAC_seppop <- seppop(QUAC_garden_wild_gen)

##calculate number of individuals per pop
n_ind_p_pop <- as.numeric(table(QUAC_seppop[[2]]@pop))

##convert the wild genind object to a genpop object
QUAC_wild_genpop <- genind2genpop(QUAC_seppop[[2]])

##create documents for comparison 
n_ind_W<-table(QUAC_garden_wild_gen@pop)[2];  n_ind_G<-table(QUAC_garden_wild_gen@pop)[1]; 
QUAC_alleles_cap <- colSums(QUAC_seppop[[1]]@tab,na.rm=T)

#######create table for % alleles captured by frequency and how many duplicates were present  
#create list with duplicates 
dup_reps <- c(0:9)

##create a table to store % alleles captured by gardens pops where no alleles are dropped 
QUAC_allele_cap_table_ndrop0 <- matrix(nrow = length(dup_reps), ncol = length(list_allele_cat))

##create a table to store % alleles captured by garden pops where alleles are dropped if there are fewer than 2
QUAC_allele_cap_table_ndrop2 <- matrix(nrow = length(dup_reps), ncol = length(list_allele_cat))

######create arrays and lists to store results 
QUAC_allele_cat <- list()
#create allele existing df
QUAC_all_exist_df <- matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))
#create df of wild alleles captured by gardens
QUAC_wild_cap_df <- matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))
##data frame to record allele capture code
QUAC_allele_cap <-matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))

##run loop to generate allelic capture table 
#the outer loop is calculating how many copies of each allele in each category exists
#the inner loop is calculating the percent capture of each allele in each frequency category 
for(ndrop in c(0,2)){
  for(dup in 1:length(dup_reps)){
    for(cat in 1:length(list_allele_cat)){
      
      if(ndrop == 0) n_drop_file <- "_dr_0"
      if(ndrop == 2) n_drop_file <- ""
      
      ##first calculate the frequency categories of alleles in the wild individuals   	
      QUAC_allele_cat <- get.allele.cat(QUAC_wild_genpop, 1, 1, as.numeric(n_ind_p_pop), n_drop = ndrop, glob_only = TRUE)	
      
      ##create a data frame with all of the alleles existing by category
      QUAC_all_exist_df[dup, cat] <- round(sum(QUAC_alleles_cap[QUAC_allele_cat[[cat]]] > dup_reps[[dup]]))
      
      ##now determine how many wild alleles were captured per category 
      QUAC_wild_cap_df[dup, cat] <- round(sum(QUAC_alleles_cap[QUAC_allele_cat[[cat]]] > dup_reps[[dup]])/length(QUAC_allele_cat[[cat]]),4)
      
      ##code to store as one data frame 
      QUAC_allele_cap[dup, cat] <- paste0(signif((QUAC_wild_cap_df[dup,cat]*100),3), "% (", signif(QUAC_all_exist_df[dup,cat],3), ")")
      
    }
  }
  
  ##format tables
  #alleles existing
  rownames(QUAC_all_exist_df) <- paste0(c(1:10), " or more copies")
  colnames(QUAC_all_exist_df) <- list_allele_cat
  #percent capture of allele types by gardens
  rownames(QUAC_wild_cap_df) <- paste0(c(1:10), " or more copies")
  colnames(QUAC_wild_cap_df) <- list_allele_cat
  #comparison of percent of wild alleles captured in garden 
  rownames(QUAC_allele_cap) <- paste0(c(1:10), " or more copies")
  colnames(QUAC_allele_cap) <- list_allele_cat

  ##write out data frames
  write.csv(QUAC_all_exist_df, paste0("C:/Users/eschumacher/Documents/GitHub/QUAC_insitu_exsitu/QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_all_exist_df", n_drop_file, ".csv"))
  write.csv(QUAC_wild_cap_df, paste0("C:/Users/eschumacher/Documents/GitHub/QUAC_insitu_exsitu/QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_wild_cap_df", n_drop_file, ".csv"))
  write.csv(QUAC_allele_cap, paste0("C:/Users/eschumacher/Documents/GitHub/QUAC_insitu_exsitu/QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_all_cap_garden_df", n_drop_file, ".csv"))
}

###write session info out
sessionInfo()
