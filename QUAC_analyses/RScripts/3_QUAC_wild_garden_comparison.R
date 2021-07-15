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



