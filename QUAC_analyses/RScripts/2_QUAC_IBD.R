################This script calculates the Fst between wild Q. acerifolia populations 
###############and compares it to the distance between populations 
##############First we calculate the Fst between wild Q. acerifolia populations 
#############and then run a linear regression between distance between populations and then run a Mantel test

##########################
######## Libraries #######
##########################

library(adegenet)
library(hierfstat)
library(geosphere)

################################
########## Load files ##########
################################
setwd("../../QUAC_data_files")

##load in genepop files as genind objects
QUAC_wild_gen <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_wild_clean.gen", ncode = 3)
##load in data frame with all scores
QUAC_wild_df <- read.csv("QUAC_data_frames/Garden_Wild/QUAC_wild_clean_df.csv")
##load in data frame with summary stats
QUAC_sumstat_df <- read.csv("../QUAC_analyses/Results/Sum_Stats/QUAC_wild_gendiv_sumstat_df.csv")

###create population name list
QUAC_pop_list <- unique(QUAC_wild_df$Pop)

###rename individuals in genind objects
rownames(QUAC_wild_gen@tab) <- QUAC_wild_df$Ind
##name populations in genind object
levels(QUAC_wild_gen@pop) <- QUAC_pop_list

#################################
###### Fst Calculations #########
#################################
##create data frame for just coordinates of each population
QUAC_coords <- QUAC_sumstat_df[,c(2:3)]

##name every pop
rownames(QUAC_coords) <- QUAC_pop_list

##convert to hierfstat format object 
QUAC_hierfstat <- genind2hierfstat(QUAC_wild_gen)

##run pairwise fst code 
QUAC_fst_df <- pairwise.neifst(QUAC_hierfstat)

##calculate geographic distances between mean locations
QUAC_dist <- matrix(nrow = length(QUAC_pop_list), ncol = length(QUAC_pop_list))

for(first in 1:length(QUAC_pop_list)){
  for(second in 1:length(QUAC_pop_list)){
    QUAC_dist[first,second] <-  distm(QUAC_coords[first,], QUAC_coords[second,], fun = distGeo)/1000
  }
}

##replacce NAs with zeroes 
QUAC_dist[is.na(QUAC_dist)] <- 0
QUAC_fst_df[is.na(QUAC_fst_df)] <- 0

##create a linear regression
QUAC_fst_dist <- lm(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)])

##visualize the isolation by distance relationship with p-value
pdf(paste0(QUAC_analysis_results, "/Clustering/QUAC_Dist_Fst.pdf"))
plot(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)], pch = 17, ylim = c(0,0.13), 
     xlim = c(0,200),
     xlab = c("Distance (km)"), ylab = c("Fst"))
abline(QUAC_fst_dist)
legend('bottomleft', legend = c("R2 = -0.12","p-value = 0.865"), bty = 'n')
dev.off()
