#These are the next steps in calculating genetic diversity 
#First, we tested linkage disequilibrium, null alleles, 
#and Hardy Weinberg equilibrium.
#Next, a data frame including expected heterozygosity, allelic richness,
#number of alleles, mean longtiude and latitude by population, and 
#individual numbers. This table is included in full in the supplemental text 
#of this manuscript.

##########################
######## Libraries #######
##########################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(diveRsity)

#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

#######################################
########## Load genind files ##########
#######################################

setwd(QUAC_data_files)

##convert to a genind 
#arp2gen(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_allpop_clean.arp"))

##load in genind for QUAC
QUAC_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_allpop_clean.gen"), ncode = 3)

##load in just garden, just wild
QUAC_garden_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_garden_clean.gen"), ncode = 3)
QUAC_wild_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_wild_clean.gen"), ncode = 3)

##load in data frame 
QUAC_allpop_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\QUAC_allpop_clean.csv"))

##rename individuals in the genind file
rownames(QUAC_gen@tab) <- QUAC_allpop_df$ID

##create population name lists
QUAC_popnames <- unique(QUAC_allpop_df$Pop)

##rename pops in genind
levels(QUAC_gen@pop) <- QUAC_popnames
levels(QUAC_garden_gen@pop) <- QUAC_garden_names
levels(QUAC_wild_gen@pop) <- QUAC_wild_names

############################################################################
####### Run Genetic Diversity Checks like LD, HWE, Null Alleles  ###########
############################################################################
##calculate null alleles 
QUAC_null_all <- null.all(QUAC_gen)

##create null allele table
QUAC_null_all_df <- matrix(nrow = length(rownames(QUAC_null_all$null.allele.freq$summary1)),
                         ncol = length(colnames(QUAC_null_all$null.allele.freq$summary1)))

QUAC_null_all_df <- QUAC_null_all$null.allele.freq$summary1

####calculate HWE deviations
##bn HWE test
QUAC_hwe_pop <- seppop(QUAC_gen) %>% lapply(hw.test, B = 1000)

##create table by populations
QUAC_HWE_allpop_df <- sapply(QUAC_hwe_pop, "[", i = TRUE, j = 3)
##name columns
colnames(QUAC_HWE_allpop_df) <- QUAC_popnames

###calculate linkage disequilibrium 
QUAC_ld <- pair.ia(QUAC_gen, sample = 1000)
QUAC_ld_df <- data.frame(round(QUAC_ld,digits = 2))

##write out null allele document 
write.csv(QUAC_null_all_df, paste0(QUAC_analysis_results, "\\Sum_Stats\\QUAC_null_all_df.csv"))
write.csv(QUAC_HWE_allpop_df, paste0(QUAC_analysis_results, "\\Sum_Stats\\QUAC_HWE_df.csv"))
write.csv(QUAC_ld_df, paste0(QUAC_analysis_results, "\\Sum_Stats\\QUAC_ld_df.csv"))

######################################
############# basic stats ############
######################################

QUAC_pop_type_gen <- list.files(path = "QUAC_genind\\garden_wild", pattern = ".gen$")
QUAC_pop_type_df <- list.files(path = "QUAC_data_frames\\garden_wild", pattern = ".csv$")
QUAC_pop_df_list <- list()
for(a in 1:2) {
  QUAC_pop_df_list[[a]] <- read.csv(paste0("QUAC_data_frames\\garden_wild\\", QUAC_pop_type_df[[a]]))
}
QUAC_garden_wild <- c("garden","wild")

for(pop in 1:length(QUAC_pop_type_gen)){
  
  QUAC_temp_gen <- read.genepop(paste0("QUAC_genind\\garden_wild\\",QUAC_pop_type_gen[[pop]]), ncode = 3)
  
  ##rename rownames
  rownames(QUAC_temp_gen@tab) <- QUAC_pop_df_list[[pop]][,1]
  ##get pop names
  QUAC_pop_names <- unique(QUAC_pop_df_list[[pop]][,2])
  ##rename pops
  levels(QUAC_temp_gen) <- QUAC_pop_names
  
  ###start genetic analyses
  QUAC_temp_sum <- summary(QUAC_temp_gen)
  ##create poppr file 
  QUAC_poppr <- poppr(QUAC_temp_gen)
  #expected heterozygosity 
  QUAC_hexp <- QUAC_poppr[1:length(QUAC_pop_names), 10]
  ##allele numbers by pop 
  QUAC_nall <- QUAC_temp_sum$pop.n.all
  ##individual numbers
  QUAC_ind <- QUAC_poppr[1:length(QUAC_pop_names), 2:3]
  ##allelic richness code 
  QUAC_alleles <- QUAC_temp_sum$pop.n.all/length(QUAC_temp_gen@loc.n.all)
  QUAC_allrich <- colMeans(allelic.richness(QUAC_temp_gen)$Ar)	
  
  ##create data frame 
  QUAC_gendiv_sumstat_df <- signif(cbind(QUAC_ind, QUAC_nall, QUAC_allrich, QUAC_hexp),3)
  
  ##name columns and rows 
  rownames(QUAC_gendiv_sumstat_df) <- QUAC_pop_names
  colnames(QUAC_gendiv_sumstat_df) <- c("Number of Individuals", "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")
  
  if(pop == 2){
    
    ##Calculate mean longitude and latitude for each wild population
    #first create matrices
    QUAC_mean_lon <- matrix()
    QUAC_mean_lat <- matrix()
    
    ##identifying mean latitudes for each population
    
    for(pop2 in QUAC_pop_names){
      
      QUAC_mean_lon[pop2] <- mean(QUAC_pop_df_list[[2]][QUAC_pop_df_list[[2]]$Pop == pop2,][,3])
      QUAC_mean_lon <- matrix(QUAC_mean_lon)
      
    }
    
    for(pop2 in QUAC_pop_names){
      
      QUAC_mean_lat[pop2] <- mean(QUAC_pop_df_list[[2]][QUAC_pop_df_list[[2]]$Pop == pop2,][,4])
      
      
    }
    
    QUAC_lonlat_df <- cbind(QUAC_mean_lon, QUAC_mean_lat)
    QUAC_lonlat_df <- QUAC_lonlat_df[c(2:6),]
    
    QUAC_gendiv_sumstat_df <- cbind(QUAC_lonlat_df, QUAC_gendiv_sumstat_df)
    
  }
  
  ##write out csv 
  write.csv(QUAC_gendiv_sumstat_df, paste0(QUAC_analysis_results, "\\Sum_Stats\\QUAC_", QUAC_garden_wild[[pop]],"_gendiv_sumstat_df.csv"))
  
}



