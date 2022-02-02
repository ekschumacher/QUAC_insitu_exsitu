##This script calculates the Fst between wild Q. acerifolia 
#compares genetic differentiation with distance
#First we calculate the Fst between wild Q. acerifolia populations 
#and then regress it with the distance between populations. 
#We use "cleaned" data files here, which means genetic 
#data files have been cleaned for clones and individuals with 
#too much missing data (25% missing data)

#########################
#        Libraries      #
#########################

library(adegenet)
library(hierfstat)
library(geosphere)

#########################
#   Load Data Files     #
#########################
#set working directory
setwd("../../QUAC_data_files")

#load in genepop files as genind objects
QUAC_wild_gen <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_wild_allpop_clean.gen", ncode = 3)

#load in data frame with summary stats
QUAC_sumstat_df <- read.csv("../QUAC_analyses/Results/Sum_Stats/Wild_gendiv_sumstat_df.csv")

#create population name list
QUAC_wildpop_names <- QUAC_sumstat_df[,1]

#name populations in genind object
levels(QUAC_wild_gen@pop) <- QUAC_wildpop_names

#############################
#     Fst Calculations      #
#############################
#create data frame for just coordinates of each population
QUAC_coords <- QUAC_sumstat_df[,c(2:3)]

#name every pop
rownames(QUAC_coords) <- QUAC_wildpop_names

#convert to hierfstat format object 
QUAC_hierfstat <- genind2hierfstat(QUAC_wild_gen)

##run pairwise fst code 
QUAC_fst_df <- pairwise.neifst(QUAC_hierfstat)

##calculate geographic distances between mean locations
QUAC_dist <- matrix(nrow = length(QUAC_wildpop_names), ncol = length(QUAC_wildpop_names))

for(first in 1:length(QUAC_wildpop_names)){
  for(second in 1:length(QUAC_wildpop_names)){
    QUAC_dist[first,second] <-  distm(QUAC_coords[first,], QUAC_coords[second,], fun = distGeo)/1000
  }
}

##replacce NAs with zeroes 
QUAC_dist[is.na(QUAC_dist)] <- 0
QUAC_fst_df[is.na(QUAC_fst_df)] <- 0

##create a linear regression
QUAC_fst_dist <- lm(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)])

##visualize the isolation by distance relationship with p-value
pdf("../QUAC_analyses/Results/Clustering/QUAC_Dist_Fst.pdf")
plot(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)], pch = 17, ylim = c(0,0.13), 
     xlim = c(0,200),
     xlab = c("Distance (km)"), ylab = c("Fst"))
abline(QUAC_fst_dist)
legend('bottomleft', legend = c("R2 = -0.12","p-value = 0.865"), bty = 'n')
dev.off()
