##This script was written to compare wild and garden diversity in 
#genomic microsatellites and possibly-gene-linked microsatellites (ESTs)
#Wild and garden individuals had their loci separated by genomic or EST 
#and then were compared for wild and garden diversity differences.
#We use "clean" data files here - these data files have had 
#clones removed and individuals with too much missing data (25%).

#####################
#     Libraries     #
#####################
library(adegenet)
library(poppr)
library(hierfstat)

###############################
#     Loading data files      #
###############################
setwd("../../QUAC_data_files")

##load in data frames for EST/genomic microsat diversity comparison
#load names of files
QUAC_EST_genepop_name_list <- list.files(path = "QUAC_adegenet_files/EST",
                                    pattern = ".gen$")

QUAC_EST_df_name_list <- list.files(path = "QUAC_data_frames/EST",
                               pattern = ".csv$")

###########################################################
#     Genomic/EST microsatellite diversity comparison     #
###########################################################
##create lists to store data files for genepop and data frame files 
QUAC_EST_genepop_list <- list()
QUAC_EST_df_list <- list()
hexp_list <- list()
allrich_list <- list()
##matrix name list 
QUAC_EST_mean_df <- matrix(nrow = 5, ncol = 2)

#pop type lists 
pop_type_list <- c("Garden_EST","Garden_nonEST", "Wild_EST", "Wild_nonEST")

##loop to write in files
for(e in 1:length(QUAC_EST_genepop_name_list)){
  
  QUAC_EST_genepop_list[[e]] <- read.genepop(paste0("QUAC_adegenet_files/EST/", 
                                                    QUAC_EST_genepop_name_list[[e]]), ncode = 3)
  ##load data frames 
  QUAC_EST_df_list[[e]] <- read.csv(paste0("QUAC_data_frames/EST/", QUAC_EST_df_name_list[[e]] ))
  
  ##name individuals within each genind object 
  rownames(QUAC_EST_genepop_list[[e]]@tab) <- QUAC_EST_df_list[[e]][,1]
  
  ##genetic diversity 
  #hexp
  hexp_list[[e]] <- data.frame(as.numeric(summary(QUAC_EST_genepop_list[[e]])[7]$Hexp))
  
  #add pop types
  hexp_list[[e]]$pop_type <- pop_type_list[[e]]
  
  #all rich
  allrich_list[[e]] <- allelic.richness(QUAC_EST_genepop_list[[e]])$Ar
  
  #add pop types 
  allrich_list[[e]]$pop_type <- pop_type_list[[e]]

  
}
#combine lists into a data frame - allelic richness 
QUAC_allrich_EST_df <- rbind(allrich_list[[1]][,c(1,3)], allrich_list[[2]][,c(1,3)], allrich_list[[3]][,c(1,3)],
                                 allrich_list[[4]][,c(1,3)])

#save mean values for allelic richness 
for(pop_type in 1:length(pop_type_list)) QUAC_EST_mean_df[pop_type,1] <- mean(QUAC_allrich_EST_df[QUAC_allrich_EST_df$pop_type == paste0(pop_type_list[[pop_type]]),][,1])

#save p-value for allelic richness 
QUAC_EST_mean_df[5,1] <- kruskal.test(QUAC_allrich_EST_df[,1]~as.factor(QUAC_allrich_EST_df[,2]))[3]$p.value
  
#combine lists into a hexp data frame 
QUAC_hexp_EST_df <- rbind(hexp_list[[1]], hexp_list[[2]], hexp_list[[3]], hexp_list[[4]])

#save mean hexp values in a data frame 
for(pop_type in 1:length(pop_type_list)) QUAC_EST_mean_df[pop_type,2] <- mean(QUAC_hexp_EST_df[QUAC_hexp_EST_df$pop_type == paste0(pop_type_list[[pop_type]]),][,1])

#save p-value
QUAC_EST_mean_df[5,2] <- kruskal.test(QUAC_hexp_EST_df[,1]~as.factor(QUAC_hexp_EST_df[,2]))[3]$p.value

#name rows and columns 
rownames(QUAC_EST_mean_df) <- c(pop_type_list, "p_value")
colnames(QUAC_EST_mean_df) <- c("All_Rich", "Hexp")
#round data frame 
QUAC_EST_mean_df <- signif(QUAC_EST_mean_df, 3)

#write out matrix 
write.csv(QUAC_EST_mean_df, "../QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_EST_mean_df.csv")

#####visualize results 
##plot boxplots
pdf("../QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_EST_boxplots.pdf")
boxplot(hexp_list[[1]][,1], hexp_list[[2]][,1], hexp_list[[3]][,1], hexp_list[[4]][,1], ylim = c(0,1),
        names = c("Garden EST","Garden Genomic", "Wild EST", "Wild Genomic"),
        main = "Expected Heterozygosity Compared Over Garden, Wild, ESTs, and Genomic Microsatellite Regions")

boxplot(allrich_list[[1]][,1], allrich_list[[2]][,1], allrich_list[[3]][,1], allrich_list[[4]][,1], 
        ylim = c(0,25), names = c("Garden EST","Garden Genomic", "Wild EST", "Wild Genomic"),
        main = "Allelic Richness Compared Over Garden, Wild, ESTs, and Genomic Microsatellite Regions")
dev.off()
