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
QUAC_garden_wild_df <- read.csv("QUAC_data_frames/QUAC_garden_wild_clean_df.csv")

##rename individuals in genind object
rownames(QUAC_garden_wild_gen@tab) <- QUAC_garden_wild_df$Ind

##rename population names in the genind object
#create pop name list
QUAC_popnames <- unique(QUAC_garden_wild_df$Pop)
#rename populations in the genind object
levels(QUAC_garden_wild_gen@pop) <- QUAC_popnames

##load in function to calculate allele frequency categories
source("../QUAC_analyses/RScripts/Fa_sample_funcs.R")

##create functions to run code 
colMax <- function(data) sapply(data, max, na.rm = TRUE)
sample.pop<-function(genind_obj,vect_pop_ID,vect_samp_sizes){
  p<-length(vect_pop_ID)
  if (p>1) {
    for (p in 1:length(vect_pop_ID))
      alleles[p,]<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)
    alleles<-colSums(alleles)
  } else {alleles<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)}
  
  alleles
}

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

##create matrices to store allele frequency/allele capture code 
QUAC_all_exist <- matrix(nrow = 1, ncol = 9)
QUAC_wild_capt <- matrix(nrow = 1, ncol = 10)

##seppop genind 
QUAC_seppop <- seppop(QUAC_garden_wild_gen)

##calculate number of individuals per pop
n_ind_p_pop <- as.numeric(table(QUAC_seppop[[2]]@pop))

##convert to a genpop object
QUAC_wild_genpop <- genind2genpop(QUAC_seppop[[2]])

##calculate 	
QUAC_allele_cat <- get.allele.cat(QUAC_wild_genpop, 1, 1, as.numeric(n_ind_p_pop), glob_only = TRUE)	

##create alleles captured by gardens 
n_ind_W<-table(QUAC_garden_wild_gen@pop)[2];  n_ind_G<-table(QUAC_garden_wild_gen@pop)[1]; 
QUAC_alleles_cap <- colSums(QUAC_seppop[[1]]@tab,na.rm=T)

##create a data frame with all of the alleles existing by category
for (i in 1:9) QUAC_all_exist[,i]<- sum((QUAC_allele_cat[[i]])> 0, na.rm=T)

##now determine how many wild alleles were captured per category 
for (j in 1:length(QUAC_allele_cat)) QUAC_wild_capt[,j]<-round(sum(QUAC_alleles_cap[QUAC_allele_cat[[j]]] > 0)/length(QUAC_allele_cat[[j]]),4)
  
##combine this into one table 
QUAC_allele_cap_table <- matrix(nrow = 1, ncol = length(list_allele_cat))

##add rownames and colnames 
colnames(QUAC_allele_cap_table) <- list_allele_cat

##run loop to generate allelic capture table 
for(m in 1:length(QUAC_all_exist[1,])){
    
  QUAC_allele_cap_table[,m] <- paste0(signif((QUAC_wild_capt[,m]*100),3), "%", "(", signif(QUAC_all_exist[,m],3), ")")

}

write.csv(QUAC_allele_cap_table, "QUAC_allele_cap_df.csv")

##################################
### Replicate Allele Capture #####
##################################
##create a table to determine allelic capture in each category with duplication 
#create list with duplicates 
dup_reps <- c(0:10)
#create a data frame for alleles existing
QUAC_all_exist_dup <- matrix(nrow = length(dup_reps), ncol = length(QUAC_allele_cat))
##create a duplicate data frame captured 
QUAC_wild_cap_dup <- matrix(nrow = length(dup_reps), ncol = length(QUAC_allele_cat))
##create a data frame with combo of the two 
QUAC_wild_cap_all_dup_df <- matrix(nrow = length(dup_reps), ncol = length(QUAC_allele_cat))
#run this code through a loop code 
for (j in 1:length(QUAC_allele_cat)) {
  for(k in 1:length(dup_reps)){
    
    ##first calculate how many alleles exist in each category in wild pops 
    QUAC_all_exist_dup[k,j]<- sum(QUAC_alleles_cap[QUAC_allele_cat[[j]]] > dup_reps[[k]])
    ##second, calculate the percentages of duplicate alleles captured in gardens 
    QUAC_wild_cap_dup[k,j] <- sum(QUAC_alleles_cap[QUAC_allele_cat[[j]]] > dup_reps[[k]])/length(QUAC_allele_cat[[j]])
    ##finally, combine dfs
    QUAC_wild_cap_all_dup_df[k,j] <- paste0(signif(QUAC_wild_cap_dup[k,j]*100, 3), "% (", signif(QUAC_all_exist_dup[k,j], 3), ")")
  }
}
##update table 
colnames(QUAC_wild_cap_all_dup_df) <- list_allele_cat
rownames(QUAC_wild_cap_all_dup_df) <- c("One Copy", "> 1 Copy", "> 2 Copies", "> 3 Copies", 
                                        "> 4 Copies", "> 5 Copies", "> 6 Copies", "> 7 Copies", "> 8 Copies",
                                        "> 9 Copies", "> 10 Copies")

##write out duplicate table 
write.csv(QUAC_wild_cap_all_dup_df, "QUAC_wild_cap_all_dup_df.csv")

##write out session info
sessionInfo()
