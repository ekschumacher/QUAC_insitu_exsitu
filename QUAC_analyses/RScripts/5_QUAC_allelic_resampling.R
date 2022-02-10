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

#######################
#     Load Files      #
#######################
setwd("../../QUAC_data_files")

#load in wild genind 
QUAC_garden_wild_genind <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_garden_wild_clean.gen", ncode = 3)

#get functions
colMax <- function(data) sapply(data, max, na.rm = TRUE)
source("../QUAC_analyses/RScripts/Fa_sample_funcs.R")

############################
#     Resampling code      #
############################
#limit genind files to just wild pop 
QUAC_wild_genind <- seppop(QUAC_garden_wild_genind)[[2]]

#set number of resampling reps 
num_reps<-1000

##list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

##now create summary arrays for n to drop = 0 and n to drop = 2
summ_results_tree_ndrop2 <-array(dim=c(nrow(QUAC_wild_genind@tab)-1, 9, num_reps)) 
summ_results_tree_ndrop0 <- array(dim=c(nrow(QUAC_wild_genind@tab)-1, 9, num_reps)) 
#list of total number of individuals
n_total_indivs <- length(QUAC_wild_genind@tab[,1])
#list of all populations
n_ind_p_pop<-table(QUAC_wild_genind@pop)
#calculate total allele frequencies
allele_freqs <- colSums(QUAC_wild_genind@tab)/(n_total_indivs*2) 
#create list for inclusion of very rare alleles or not
n_drop <- c(0,2)

#loop to resample wild trees and determine at what sample size alleles are calculated 
#outer loop runs it for inclusion of very rare alleles or not (n to drop = 0 includes all rare alleles)
#n to drop = 2 removes very rare alleles 
#inner loop calculates sample size for allelic resampling 
for(ndrop in n_drop){
  #Repeat the resampling many times
  for (nrep in 1:num_reps){ 
  ##calculate the allelic category 
  allele_cat <- get.allele.cat(QUAC_wild_genind, region_makeup=NULL, 2, n_ind_p_pop, n_drop=ndrop, glob_only=T)

  #create empty matrix
  alleles_samp <- matrix(nrow=nrow(QUAC_wild_genind@tab)-1,ncol=length(allele_freqs))
  
  #This loop will sample trees from t = 2 to the total number of trees
  for (t in 2:(nrow(QUAC_wild_genind@tab)-1)){ ##minus one because our number of trees is 172, have to compare between 2 trees 
    
    #create a sample of trees of length t, by using 'sample()' which randomly samples rows
    alleles_samp <- colSums(QUAC_wild_genind@tab[sample(1:nrow(QUAC_wild_genind@tab), t),],na.rm=T)
    
    ##now add a loop to store results 
    if(ndrop == 2){
      
      #store allele frequency summaries in an array - for dropping very rare alleles
      for (l in 1:length(allele_cat)) summ_results_tree_ndrop2[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
      
    }else{
      #store allele frequency summaries in an array - including very rare alleles 
      for (l in 1:length(allele_cat)) summ_results_tree_ndrop0[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
    }
   }
  }  
}    

##Divide by numbers of alleles in each category to calculate % of each frequency of alleles captured
#for n to drop = 2
for (n in 1:num_reps) summ_results_tree_ndrop2[,,n]<-t(t(summ_results_tree_ndrop2[,,n])/summ_results_tree_ndrop2[length(summ_results_tree_ndrop2[,1,1]),,n])
#for n to drop = 0 
for (n in 1:num_reps) summ_results_tree_ndrop0[,,n]<-t(t(summ_results_tree_ndrop0[,,n])/summ_results_tree_ndrop0[length(summ_results_tree_ndrop0[,1,1]),,n])

########################################
#     Reporting Resampling Results     #
########################################
#create a list to store frequency summary results
all_mean_list <- list()

#create a data frame with mean  
all_mean_list[[1]] <- as.data.frame(apply(summ_results_tree_ndrop0[,,1:num_reps],c(1,2),mean,na.rm=T)*100)[-1,]
all_mean_list[[2]] <- as.data.frame(apply(summ_results_tree_ndrop2[,,1:num_reps],c(1,2),mean,na.rm=T)*100)[-1,]

##loop to plot resampling results including and removing very rare alleles  
for(ndrop in 1:length(all_mean_list)){
  
  #write PDF with name
  pdf(paste0("../QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_all_resampling_ndrop", n_drop[[ndrop]],".pdf"))
  #add points
  plot(all_mean_list[[ndrop]][,1], col = "red", pch = 20, xlab = "Number of Individuals", 
       ylab = "Percent Diversity Capture", xlim = c(0,171), ylim = c(0,100), cex = 1.2,
       main = "Percent Diversity Capture (All Alleles Included)")
  points(all_mean_list[[ndrop]][,2], col = "firebrick", pch = 20, cex = 1.2)
  points(all_mean_list[[ndrop]][,3], col = "darkorange3", pch = 20, cex = 1.2)
  points(all_mean_list[[ndrop]][,4], col = "coral", pch = 20, cex = 1.2)
  points(all_mean_list[[ndrop]][,5], col = "deeppink4", pch = 20, cex = 1.2)
  
  dev.off()
}

##Create data frame of min sample size to sample 95% of diversity
#create a data frame to store results 
min_samp_95 <- matrix(nrow = length(n_drop), ncol = length(list_allele_cat))

##loop to calculate min sample size
for(ndrop in 1:length(n_drop)){ 
 for(col in 1:length(list_allele_cat)){
  
    min_samp_95[ndrop,col] <- which(all_mean_list[[ndrop]][,col] >= 95)[1]
  
 }
}
##name rows and columns of the matrix 
rownames(min_samp_95) <- c("N Ind for 95% (All Alleles)", 
                           "N Ind for 95% (Rare Alleles Dropped)")
colnames(min_samp_95) <- list_allele_cat

##write out data frame 
write.csv(min_samp_95, "../QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_min_samp_95.csv")
