#############This script was written to compare wild and garden diversity in 
###########genomic microsatellites and possibly-gene-linked microsatllites (ESTs)
#########Wild and garden individuals had their loci separated by genomic or EST 
########and then were compared for wild and garden diversity differences

#######################
###### Libraries ######
#######################
library(adegenet)
library(poppr)
library(hierstat)

######################################
### Loop for comparing diversity - wild/garden, EST vs. genomic
##############################
setwd("../../QUAC_data_files")

####load in data frames for EST/genomic microsat diversity comparison
#load names of files
QUAC_EST_genepop_name_list <- list.files(path = "QUAC_adegenet_files/EST",
                                    pattern = ".gen$")

QUAC_EST_df_name_list <- list.files(path = "QUAC_data_frames/EST",
                               pattern = ".csv$")

##create lists to store data files for genepop and data frame files 
QUAC_EST_genepop_list <- list()
QUAC_EST_df_list <- list()
hexp_list <- list()
allrich_list <- list()
##matrix name list 
QUAC_EST_mean_df <- matrix(nrow = 5, ncol = 2)

##loop to write in files
for(e in 1:length(QUAC_EST_list)){
  
  QUAC_EST_genepop_list[[e]] <- read.genepop(paste0("QUAC_adegenet_files/EST/", 
                                                    QUAC_EST_list[[e]]), ncode = 3)
  ##load data frames 
  QUAC_EST_df_list[[e]] <- read.csv(paste0("QUAC_data_frames/EST/", QUAC_EST_df_name_list[[e]] ))
  
  ##name individuals within each genind object 
  rownames(QUAC_EST_genepop_list[[e]]@tab) <- QUAC_EST_df_list[[e]][,1]
  
  hexp_list[[e]] <- as.numeric(summary(QUAC_EST_genepop_list[[e]])[7]$Hexp)
  
  allrich_list[[e]] <- allelic.richness(QUAC_EST_genepop_list[[e]])$Ar
  
  ##generate means table 
  #add allrich calculations
  QUAC_EST_mean_df[e,1] <- mean(allrich_list[[e]][,1])
  #add hexp meeans
  QUAC_EST_mean_df[e,2] <- mean(hexp_list[[e]])
  
  
}

##create matrix to store results 
QUAC_EST_gendiv_comp <- matrix(nrow = 30, ncol = 3)
##paste in individuals 
QUAC_EST_gendiv_comp[1:10,1] <- "QUAC_garden_EST"
QUAC_EST_gendiv_comp[11:15,1] <- "QUAC_garden_non_EST"
QUAC_EST_gendiv_comp[16:25,1] <- "QUAC_wild_EST"
QUAC_EST_gendiv_comp[26:30,1] <- "QUAC_wild_non_EST"
##now paste in values 
QUAC_EST_gendiv_comp[1:10,2] <- as.numeric(allrich_list[[1]][,1])
QUAC_EST_gendiv_comp[11:15,2] <- allrich_list[[2]][,1]
QUAC_EST_gendiv_comp[16:25,2] <- allrich_list[[3]][,1]
QUAC_EST_gendiv_comp[26:30,2] <- allrich_list[[4]][,1]
##paste in hexp values 
QUAC_EST_gendiv_comp[1:10,3] <- hexp_list[[1]]
QUAC_EST_gendiv_comp[11:15,3] <- hexp_list[[2]]
QUAC_EST_gendiv_comp[16:25,3] <- hexp_list[[3]]
QUAC_EST_gendiv_comp[26:30,3] <- hexp_list[[4]]

##save p-value 
QUAC_EST_allrich_pvalue <- kruskal.test(as.numeric(QUAC_EST_gendiv_comp[,2])~as.factor(QUAC_EST_gendiv_comp[,1]))[3]
QUAC_EST_hexp_pvalue <- kruskal.test(as.numeric(QUAC_EST_gendiv_comp[,3])~as.factor(QUAC_EST_gendiv_comp[,1]))[3]

##add p-value to the mean table 
QUAC_EST_mean_df[5,1] <- unlist(QUAC_EST_allrich_pvalue$p.value)
QUAC_EST_mean_df[5,2] <- unlist(QUAC_EST_hexp_pvalue$p.value)
##now add information and write out matrix 
rownames(QUAC_EST_mean_df) <- c(sub("\\..*","",QUAC_EST_df_name_list), "p-value")
colnames(QUAC_EST_mean_df) <- c("Allelic Richness", "Expected Heterozygosity")
##round
QUAC_EST_mean_df <- signif(QUAC_EST_mean_df, 3)
##write out matrix 
write.csv(QUAC_EST_mean_df, "../QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_EST_mean_df.csv")
