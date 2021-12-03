##############This script details the analysis of actually comparing genetic diversity levels between garden 
#############and wild Q. acerifolia populations. We first calculated diversity levels throughout all garden  
############and wild populations (indicated by allelic richness and expected heterozygosity) and then ran a t-test 
###########between the values. The resulting figures can be found in the main text 

##########################
######## Libraries #######
##########################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)

################################
########## Load files ##########
################################
setwd("../../QUAC_data_files")

##convert to a genepop file, if necessary
#arp2gen(paste0(QUAC_data_files, "/QUAC_genind/QUAC_garden_wild_clean.arp"))

##load in genepop file as a genind object 
QUAC_garden_wild_gen <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_garden_wild_clean.gen", ncode = 3)

##load in data frame 
QUAC_garden_wild_df <- read.csv("QUAC_data_frames/Garden_Wild/QUAC_garden_wild_clean_df.csv")

##rename individuals in genind object
rownames(QUAC_garden_wild_gen@tab) <- QUAC_garden_wild_df$ID

##rename population names in the genind object
#create pop name list
QUAC_popnames <- unique(QUAC_garden_wild_df$POP)
#rename populations in the genind object
levels(QUAC_garden_wild_gen@pop) <- QUAC_popnames

##load in function to calculate allele frequency categories
source("../QUAC_analyses/RScripts/Fa_sample_funcs.R")

##create functions to run code 
colMax <- function(data) sapply(data, max, na.rm = TRUE)

############################################################
############ Comparing wild and garden populations #########
############################################################
setwd("../QUAC_analyses/Results/Wild_Garden_Comparison")
##calculate genetic diversity statistics 
QUAC_wild_hexp <- summary(seppop(QUAC_garden_wild_gen)$Wild)[7]
QUAC_garden_hexp <- summary(seppop(QUAC_garden_wild_gen)$Garden)[7]

##calculate allelic richness
QUAC_wild_allrich <- allelic.richness(seppop(QUAC_garden_wild_gen)$Wild)$Ar
QUAC_garden_allrich <- allelic.richness(seppop(QUAC_garden_wild_gen)$Garden)$Ar

##create gendiv data frame 
QUAC_gendiv_df <- matrix(nrow = 30, ncol = 3)
QUAC_gendiv_df[1:15,1] <- c("Garden")
QUAC_gendiv_df[16:30,1] <- c("Wild")
##load in heterozygosity
QUAC_gendiv_df[1:15,2] <- QUAC_garden_hexp$Hexp
QUAC_gendiv_df[16:30,2] <- QUAC_wild_hexp$Hexp
##load in allelic richness 
QUAC_gendiv_df[1:15,3] <- QUAC_wild_allrich[,1]
QUAC_gendiv_df[16:30,3] <- QUAC_garden_allrich[,1]

#####test for normality and homogeneity of variance --hexp
##test homogeneity of variances 
var.test(as.numeric(QUAC_gendiv_df[,2])~as.factor(QUAC_gendiv_df[,1]))
##normality test
shapiro.test(as.numeric(QUAC_gendiv_df[,2]))
##run t-test 
hexp_krusk <- kruskal.test(as.numeric(QUAC_gendiv_df[,2])~as.factor(QUAC_gendiv_df[,1]))

#####t-test for allelic richness 
##test for homogeneity of variance
var.test(as.numeric(QUAC_gendiv_df[,3])~as.factor(QUAC_gendiv_df[,1]))
##test for normal data 
shapiro.test(as.numeric(QUAC_gendiv_df[,3]))
##do kruskal-wallis
##run t-test 
allrich_krusk <- kruskal.test(as.numeric(QUAC_gendiv_df[,3])~as.factor(QUAC_gendiv_df[,1]))

##results data frame 
QUAC_garden_wild_results <- matrix(nrow = 2, ncol = 3)
##name columns and rows 
rownames(QUAC_garden_wild_results) <- c("Expected Heterozygosity", "Allelic Richness")
colnames(QUAC_garden_wild_results) <- c("Garden Mean", "Wild Mean", "P-Value")

##add values 
QUAC_garden_wild_results[1,3] <- hexp_krusk[3]$p.value
QUAC_garden_wild_results[2,3] <- allrich_krusk[3]$p.value
QUAC_garden_wild_results[1,1] <- mean(QUAC_garden_hexp$Hexp) ; QUAC_garden_wild_results[1,2] <- mean(QUAC_wild_hexp$Hexp)
QUAC_garden_wild_results[2,1] <- mean(QUAC_garden_allrich[,1]) ; QUAC_garden_wild_results[2,2] <- mean(QUAC_wild_allrich[,1])

##write out table 
write.csv(QUAC_garden_wild_results, "QUAC_garden_wild_gendiv_dif_df.csv")

####################################
#### Allelic capture code ##########
####################################
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
      QUAC_all_exist_df[dup, cat] <- sum((QUAC_allele_cat[[cat]])> dup_reps[[dup]])
      
      ##now determine how many wild alleles were captured per category 
      QUAC_wild_cap_df[dup, cat] <- round(sum(QUAC_alleles_cap[QUAC_allele_cat[[cat]]] > dup_reps[[dup]])/length(QUAC_allele_cat[[cat]]),4)
      
      ##code to store as one data frame 
      #QUAC_allele_cap[dup, cat] <- paste0(signif((QUAC_wild_cap_df[dup,cat, ndrop]),3), "(", signif(QUAC_all_exist_df[dup,cat, ndrop],3), ")")
      
    }
  }
  ##format tables
  #alleles existing
  rownames(QUAC_all_exist_df) <- paste0(c(1:10), " or more copies")
  colnames(QUAC_all_exist_df) <- list_allele_cat
  #percent capture of allele types by gardens
  rownames(QUAC_wild_cap_df) <- paste0(c(1:10), " or more copies")
  colnames(QUAC_wild_cap_df) <- list_allele_cat
  
  ##write out data frames
  write.csv(QUAC_all_exist_df, paste0("QUAC_all_exist_df", n_drop_file, ".csv"))
  write.csv(QUAC_wild_cap_df, paste0("QUAC_wild_cap_df", n_drop_file, ".csv"))
}

###write session info out
sessionInfo()