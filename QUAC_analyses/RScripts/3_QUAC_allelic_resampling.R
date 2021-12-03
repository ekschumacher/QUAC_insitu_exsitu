num_reps<-100
colMax <- function(data) sapply(data, max, na.rm = TRUE)
source("Fa_sample_funcs.R")

##write in reps
num_reps<-1000
##begin resampling code 
QUAC_wild_genind<- QUAC_seppop[[2]]
##now create summary array 
summ_results_tree<-array(dim=c(nrow(QUAC_wild_genind@tab)-1, 11, num_reps)) ##why 11? 11 allele copies? 
n_total_indivs <- length(QUAC_wild_genind@tab[,1])
n_ind_p_pop<-table(QUAC_wild_genind@pop)
allele_freqs <- colSums(QUAC_wild_genind@tab)/(n_total_indivs*2) ##times 2 because it is a diploid
allele_cat <- get.allele.cat(QUAC_wild_genind, region_makeup=NULL, 2, n_ind_p_pop, n_drop=2, glob_only=T)

#Repeat the resampling many times
for (nrep in 1:num_reps) {
  #create empty matrix
  alleles_samp <- matrix(nrow=nrow(QUAC_wild_genind@tab)-1,ncol=length(allele_freqs))
  #This loop will sample trees from t = 2 to the total number of trees
  for (t in 2:(nrow(QUAC_wild_genind@tab)-1)){ ##minus one because our number of trees is 172, have to compare between 2 trees 
    #create a sample of trees of length t, by using 'sample()' which randomly samples rows
    alleles_samp <- colSums(QUAC_wild_genind@tab[sample(1:nrow(QUAC_wild_genind@tab), t),],na.rm=T)
    #Then simply compare that sample to your wild population with allele_cat
    for (l in 1:length(allele_cat)) summ_results_tree[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
  }
}

#Divide by the number of alleles
for (n in 1:num_reps) summ_results_tree[,,n]<-t(t(summ_results_tree[,,n])/summ_results_tree[length(summ_results_tree[,1,1]),,n])
#mean across reps using apply
all_mean<-apply(summ_results_tree[,,1:num_reps],c(1,2),mean,na.rm=T)*100

all_mean <- as.data.frame(all_mean)
all_mean <- all_mean[,1:9]
colnames(all_mean) <- list_allele_cat
all_mean <- all_mean[-1,]

##write out resampling code 
write.csv(all_mean, "allelic_resampling_table.csv")
pdf("all_resampling.pdf")

plot(all_mean[,1], col = "red", pch = 20, xlab = "Number of Individuals", 
     ylab = "Percent Diversity Capture", xlim = c(0,171), ylim = c(0,100), cex = 1.2)
points(all_mean[,2], col = "firebrick", pch = 20, cex = 1.2)
points(all_mean[,3], col = "darkorange3", pch = 20, cex = 1.2)
points(all_mean[,4], col = "coral", pch = 20, cex = 1.2)
points(all_mean[,5], col = "deeppink4", pch = 20, cex = 1.2)

legend('bottomright', legend = c("Global Alleles", "Globally Very Common Alleles", 
                                 "Globally Common", "Global Low Frequency Alleles",
                                 "Globally Rare Alleles"), pch = 20,
       col = c("red", "firebrick", "darkorange3", "coral", "deeppink4"), pt.cex = 2)
dev.off()
