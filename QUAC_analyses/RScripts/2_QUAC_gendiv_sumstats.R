###These are the next steps in calculating genetic diversity 
##First, we tested linkage disequilibrium, null alleles, and Hardy Weinberg equilibrium.
#Next, a data frame including expected heterozygosity, allelic richness, number of alleles, 
#mean longtitude and latitude by population, and individual numbers. 
#This table is included in full in the supplemental text of this manuscript.

#########################
#        Libraries      #
#########################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)

#########################
#   Load Data Files     #
#########################
##set working directory to load in data files
setwd("../../QUAC_data_files")

##load in QUAC all pop genepop as a genind object
QUAC_allpop_gen <- read.genepop("QUAC_adegenet_files/QUAC_allpop_clean.gen", ncode = 3)

#load all pop data frames
QUAC_allpop_df <- read.csv("QUAC_data_frames/QUAC_allpop_clean.csv")

##create separate lists for wild and garden population names 
QUAC_gardenpop_names <- QUAC_allpop_names[1:17]
QUAC_wildpop_names <- QUAC_allpop_names[18:22]

##load in lon/lat data frame for wild pops
QUAC_wild_lonlat_clean_df <- read.csv("QUAC_data_frames/QUAC_wild_lonlat_allpop_clean_df.csv")

############################################################
#  Null Alleles, HWE Deviation, Linkage Disequilibrium     #
############################################################
####calculate null alleles 
QUAC_nullall_pop <- lapply(QUAC_allpop_gen, null.all)
#create null allele table
QUAC_null_all_df <- signif(data.frame(QUAC_nullall_pop$null.allele.freq$summary2),3)

####calculate HWE deviations
##bn HWE test
QUAC_hwe_pop <- seppop(QUAC_allpop_gen) %>% lapply(hw.test, B = 1000)
##create table by populations
QUAC_HWE_allpop_df <- sapply(QUAC_hwe_pop, "[", i = TRUE, j = 3)
##name columns
colnames(QUAC_HWE_allpop_df) <- QUAC_popnames
##round to the 3rd digit
QUAC_HWE_allpop_df <- signif(QUAC_HWE_allpop_df, 3)

###calculate linkage disequilibrium 
QUAC_ld <- pair.ia(QUAC_allpop_gen, sample = 1000)
QUAC_ld_df <- data.frame(round(QUAC_ld,digits = 2))

##write out null allele document 
write.csv(QUAC_null_all_df, "../QUAC_analyses/Results/Sum_Stats/QUAC_null_all_df.csv")
write.csv(QUAC_HWE_allpop_df, "../QUAC_analyses/Results/Sum_Stats/QUAC_HWE_df.csv")
write.csv(QUAC_ld_df, "../QUAC_analyses/Results/Sum_Stats/QUAC_ld_df.csv")

###########################################
#          Genetic Stats by Pop           #
###########################################
##Garden vs. wild diversity capture 
#create lists for loop of the object and files used in analyses
QUAC_pop_type_gen <- list.files(path = "QUAC_adegenet_files/Garden_Wild", pattern = "allpop_clean.gen")
QUAC_pop_type_df <- list.files(path = "QUAC_data_frames/Garden_Wild", pattern = "allpop_clean_df.csv")
#create list of 'pop type' for naming output files
QUAC_pop_type <- c("QUAC_garden","QUAC_wild")

##create data frame list 
QUAC_pop_df_list <- list()

##first read in data frames in
for(a in 1:length(QUAC_pop_type_df)) {
  
  QUAC_pop_df_list[[a]] <- read.csv(paste0("QUAC_data_frames/Garden_Wild/", QUAC_pop_type_df[[a]]))

}

##create a list to save all the data frames 
QUAC_allpop_gendiv_sumstat_list <- list()

##create lists for genetic diversity statistics 
QUAC_allrich_list <- list()
QUAC_hexp_list <- list()

##write loop to calculate all summary stats 
for(pop in 1:length(QUAC_pop_type_gen)){
 
  ##load in genepop file as a genind object
  QUAC_temp_gen <- read.genepop(paste0("QUAC_adegenet_files/Garden_Wild/", QUAC_pop_type_gen[[pop]]), ncode = 3)
  
  ##name pops  
  if(pop == 1){
    levels(QUAC_temp_gen@pop) <- QUAC_gardenpop_names
  }else{
    levels(QUAC_temp_gen@pop) <- QUAC_wildpop_names
  }
  
  ##rename rownames
  rownames(QUAC_temp_gen@tab) <- QUAC_pop_df_list[[pop]][,1]
  
  ###start genetic analyses
  QUAC_temp_sum <- summary(QUAC_temp_gen)
  ##create poppr file 
  QUAC_poppr <- poppr(QUAC_temp_gen)
  ##save mean for final output table 
  QUAC_hexp_mean <- QUAC_poppr[1:length(levels(QUAC_temp_gen@pop)),10]
  ##allele numbers by pop 
  QUAC_nall <- QUAC_temp_sum$pop.n.all
  ##individual numbers
  QUAC_ind <- QUAC_poppr[1:length(levels(QUAC_temp_gen@pop)), 2:3]
  ##save allelic richness for comparison
  QUAC_allrich_list[[pop]] <- allelic.richness(QUAC_temp_gen)$Ar
  QUAC_allrich_mean <- colMeans(allelic.richness(QUAC_temp_gen)$Ar)	
  
  ##create data frame 
  QUAC_allpop_gendiv_sumstat_df <- signif(cbind(QUAC_ind, QUAC_nall, QUAC_allrich_mean, QUAC_hexp_mean),3)
   
  ##create extra columns for mean lon/lat 
  if(pop == 2){
    
    ##Calculate mean longitude and latitude for each wild population
    #first create matrices
          QUAC_mean_lon <- matrix()
          QUAC_mean_lat <- matrix()
    
    ##identifying mean longitude and latitude for each population 
          for(pop2 in QUAC_wildpop_names){
      
              QUAC_mean_lon[pop2] <- mean(QUAC_pop_df_list[[pop]][QUAC_pop_df_list[[pop]]$Pop == pop2,][,3])
              QUAC_mean_lon <- matrix(QUAC_mean_lon)
      
       }
          for(pop2 in QUAC_wildpop_names){
      
                QUAC_mean_lat[pop2] <- mean(QUAC_pop_df_list[[pop]][QUAC_pop_df_list[[pop]]$Pop == pop2,][,4])
               QUAC_mean_lat <- matrix(QUAC_mean_lat)
      
       }
    
      QUAC_lonlat_df <- cbind(QUAC_mean_lon, QUAC_mean_lat)
      QUAC_lonlat_df <- QUAC_lonlat_df[c(2:6),]
    
      QUAC_allpop_gendiv_sumstat_df <- cbind(QUAC_lonlat_df, QUAC_allpop_gendiv_sumstat_df)
      
      ##name columns 
      colnames(QUAC_allpop_gendiv_sumstat_df) <- c("Mean_Lon", "Mean_Lat", "N", "MLG", "Number of Alleles", "Allelic Richness", "Expected Heterozygosity")
      rownames(QUAC_allpop_gendiv_sumstat_df) <- QUAC_wildpop_names
    
  }else{
    ##name columns 
    colnames(QUAC_allpop_gendiv_sumstat_df) <- c("N", "MLG", "Number of Alleles", "Allelic Richness", "Expected Heterozygosity")
    rownames(QUAC_allpop_gendiv_sumstat_df) <- QUAC_gardenpop_names
     }
  
  ##write out csv 
 write.csv(QUAC_allpop_gendiv_sumstat_df, paste0("../QUAC_analyses/Results/Sum_Stats/", QUAC_pop_type[[pop]], "_gendiv_sumstat_df.csv"))
  
}
