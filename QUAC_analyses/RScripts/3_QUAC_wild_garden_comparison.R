##########################
######## Libraries #######
##########################

library(adegenet)
library(diveRsity)
library(poppr)

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

##load in get allele category function
source("C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\RScripts\\Fa_sample_funcs.R")

##load in functions
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
##calculate genetic diversity statistics 
QUAC_wild_hexp <- summary(seppop(QUAC_wild_garden_gen)$Wild)[7]
QUAC_garden_hexp <- summary(seppop(QUAC_wild_garden_gen)$Garden)[7]

##calculate allelic richness
QUAC_wild_allrich <- allelic.richness(seppop(QUAC_wild_garden_gen)$Wild)$Ar
QUAC_garden_allrich <- allelic.richness(seppop(QUAC_wild_garden_gen)$Garden)$Ar

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
write.csv(QUAC_garden_wild_results, paste0(QUAC_analysis_results, "\\Wild_Garden_Comparison\\gendiv_dif.csv"))

####################################
#### Allelic capture code ##########
####################################
##list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

##create matrices to store results
QUAC_all_exist <- matrix(nrow = 1, ncol = 9)
QUAC_wild_capt <- matrix(nrow = 1, ncol = 10)

##seppop genind 
QUAC_seppop <- seppop(QUAC_wild_garden_gen)

##calculate number of individuals per pop
n_ind_p_pop <- table(QUAC_wild_garden_gen@pop)

##convert to genpop
QUAC_wild_garden_genpop <- genind2genpop(QUAC_wild_garden_gen)

# allele_cat<-get.allele.cat(Spp_genpop, region_makeup, 2, n_ind_p_pop,n_drop=n_to_drop)	
QUAC_allele_cat <- get.allele.cat(QUAC_wild_garden_genpop, c(1:2), 2, n_ind_p_pop)	

##create alleles captured by gardens 
n_ind_W<-table(QUAC_wild_garden_gen@pop)[2];  n_ind_G<-table(QUAC_wild_garden_gen@pop)[1]; 
QUAC_alleles_cap <- colSums(QUAC_seppop[[1]]@tab,na.rm=T)

##create a data frame with all of the alleles existing by species
for (i in 1:9) QUAC_all_exist[,i]<- (sum((QUAC_allele_cat[[i]])>0,na.rm=T))

##now determine how many wild alleles were captured per category 
for (j in 1:length(QUAC_allele_cat)) QUAC_wild_capt[,j]<-round(sum(QUAC_alleles_cap[QUAC_allele_cat[[j]]]>0)/length(QUAC_allele_cat[[j]]),4)

##add to this data frame the number of individuals in gardens
QUAC_wild_capt[,10] <- n_ind_G

##combine this into one table 
QUAC_allele_cap_table <- matrix(nrow = 1, ncol = length(QUAC_all_exist[1,]))

##add rownames and colnames 
colnames(QUAC_allele_cap_table) <- list_allele_cat

##run loop to generate table 
for(m in 1:length(QUAC_all_exist[1,])){
    
  QUAC_allele_cap_table[,m] <- paste0(signif((QUAC_wild_capt[,m]*100),3), "%", "(", signif(QUAC_all_exist[,m],3), ")")

}

write.csv(QUAC_allele_cap_table, paste0(QUAC_analysis_results, "\\Wild_Garden_Comparison\\QUAC_allele_cap_table.csv"))

