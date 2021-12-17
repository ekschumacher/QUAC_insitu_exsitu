#############################
###### Libraries ############
#############################
library(adegenet)

#############################
###### Load Files ###########
setwd("../../QUAC_data_files")

##load in wild genind 
QUAC_wild_genind <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_wild_clean.gen", ncode = 3)

##get functions
colMax <- function(data) sapply(data, max, na.rm = TRUE)
source("../QUAC_analyses/RScripts/Fa_sample_funcs.R")

#############################
###### Start rep code #######
#############################
##set number of resampling reps 
num_reps<-1000

##list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

##now create summary array 
summ_results_tree_ndrop2 <-array(dim=c(nrow(QUAC_wild_genind@tab)-1, 9, num_reps)) 
summ_results_tree_ndrop0 <- array(dim=c(nrow(QUAC_wild_genind@tab)-1, 9, num_reps)) 
n_total_indivs <- length(QUAC_wild_genind@tab[,1])
n_ind_p_pop<-table(QUAC_wild_genind@pop)
allele_freqs <- colSums(QUAC_wild_genind@tab)/(n_total_indivs*2) 

n_drop <- c(0,2)

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
      
      for (l in 1:length(allele_cat)) summ_results_tree_ndrop2[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
      
      #Divide by the number of alleles
      for (n in 1:num_reps) summ_results_tree_ndrop2[,,n]<-t(t(summ_results_tree_ndrop2[,,n])/summ_results_tree_ndrop2[length(summ_results_tree_ndrop2[,1,1]),,n])
      
      
    }else{
        
        for (l in 1:length(allele_cat)) summ_results_tree_ndrop0[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
        
        for (n in 1:num_reps) summ_results_tree_ndrop0[,,n]<-t(t(summ_results_tree_ndrop0[,,n])/summ_results_tree_ndrop0[length(summ_results_tree_ndrop0[,1,1]),,n])
        
          
    }
    
   
    
   }
  }  
}    

#################################
#### Collate Sampling Means #####
#################################
###ndrop0 
#mean across reps using apply
all_mean_ndrop0 <- apply(summ_results_tree_ndrop0[,,1:num_reps],c(1,2),mean,na.rm=T)*100
all_mean_ndrop0 <- as.data.frame(all_mean_ndrop0)
all_mean_ndrop0 <- all_mean_ndrop0[,1:9]
colnames(all_mean_ndrop0) <- list_allele_cat
all_mean_ndrop0 <- all_mean_ndrop0[-1,]

##write out resampling code 
write.csv(all_mean_ndrop0, "../QUAC_analyses/Results/Wild_Garden_Comparison/allelic_resampling_table.csv")


###plot resampling graph 
pdf("../QUAC_analyses/Results/Wild_Garden_Comparison/all_resampling_ndrop0.pdf")

plot(all_mean_ndrop0[,1], col = "red", pch = 20, xlab = "Number of Individuals", 
     ylab = "Percent Diversity Capture", xlim = c(0,171), ylim = c(0,100), cex = 1.2,
     main = "Percent Diversity Capture (All Alleles Included)")
points(all_mean_ndrop0[,2], col = "firebrick", pch = 20, cex = 1.2)
points(all_mean_ndrop0[,3], col = "darkorange3", pch = 20, cex = 1.2)
points(all_mean_ndrop0[,4], col = "coral", pch = 20, cex = 1.2)
points(all_mean_ndrop0[,5], col = "deeppink4", pch = 20, cex = 1.2)

legend('bottomright', legend = c("Global Alleles", "Globally Very Common Alleles", 
                                 "Globally Common", "Global Low Frequency Alleles",
                                 "Globally Rare Alleles"), pch = 20,
       col = c("red", "firebrick", "darkorange3", "coral", "deeppink4"), pt.cex = 2)
dev.off()

##########Now plot allelic means for dropping rare alleles 
#mean across reps using apply
all_mean_ndrop2 <- apply(summ_results_tree_ndrop2[,,1:num_reps],c(1,2),mean,na.rm=T)*100
all_mean_ndrop2 <- as.data.frame(all_mean_ndrop2)
all_mean_ndrop2 <- all_mean_ndrop2[,1:9]
colnames(all_mean_ndrop2) <- list_allele_cat
all_mean_ndrop2 <- all_mean_ndrop2[-1,]

###plot resampling graph 
pdf("../QUAC_analyses/Results/Wild_Garden_Comparison/all_resampling_ndrop2.pdf")

plot(all_mean_ndrop2[,1], col = "red", pch = 20, xlab = "Number of Individuals", 
     ylab = "Percent Diversity Capture", xlim = c(0,171), ylim = c(0,100), cex = 1.2,
     main = "Percent Diversity Capture (Rare Alleles Dropped)")
points(all_mean_ndrop0[,2], col = "firebrick", pch = 20, cex = 1.2)
points(all_mean_ndrop0[,3], col = "darkorange3", pch = 20, cex = 1.2)
points(all_mean_ndrop0[,4], col = "coral", pch = 20, cex = 1.2)
points(all_mean_ndrop0[,5], col = "deeppink4", pch = 20, cex = 1.2)

legend('bottomright', legend = c("Global Alleles", "Globally Very Common Alleles", 
                                 "Globally Common", "Global Low Frequency Alleles",
                                 "Globally Rare Alleles"), pch = 20,
       col = c("red", "firebrick", "darkorange3", "coral", "deeppink4"), pt.cex = 2)
dev.off()
