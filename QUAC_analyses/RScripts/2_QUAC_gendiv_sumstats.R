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
arp2gen(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_allpop_clean.arp"))

##load in genind for QUAC
QUAC_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_allpop_clean.gen"), ncode = 3)

##load in data frame 
QUAC_allpop_clean <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\QUAC_allpop_clean.csv"))

##rename individuals in the genind file
rownames(QUAC_gen@tab) <- QUAC_allpop_clean$ID

##rename population names 
QUAC_popnames <- unique(QUAC_allpop_clean$Pop)

levels(QUAC_gen@pop) <- QUAC_popnames

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
##reorg data file 
QUAC_summary_gen <- summary(QUAC_gen)

##create poppr file 
QUAC_poppr <- poppr(QUAC_gen)

##expected heterozygosity 
QUAC_hexp <- QUAC_poppr[1:22, 10]
##allele numbers by pop 
QUAC_nall <- QUAC_summary_gen$pop.n.all
##individual numbers
QUAC_ind <- QUAC_poppr[1:22, 2:3]
##allelic richness code 
QUAC_alleles <- QUAC_summary_gen$pop.n.all/length(QUAC_gen@loc.n.all)
QUAC_allrich <- colMeans(allelic.richness(QUAC_gen)$Ar)	

##create data frame 
QUAC_gendiv_sumstat_df <- signif(cbind(QUAC_ind, QUAC_nall, QUAC_allrich, QUAC_hexp),3)

##name columns and rows 
rownames(QUAC_gendiv_sumstat_df) <- QUAC_popnames
colnames(QUAC_gendiv_sumstat_df) <- c("Number of Individuals", "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")

##write out csv 
write.csv(QUAC_gendiv_sumstat_df, paste0(QUAC_analysis_results, "\\Sum_Stats\\QUAC_gendiv_sumstat_df.csv"))
