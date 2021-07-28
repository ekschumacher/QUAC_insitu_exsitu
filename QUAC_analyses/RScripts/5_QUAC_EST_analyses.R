##This code is to compare coding regions of DNA to neutral
##And between populations

##########################
######## Libraries #######
##########################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)

#####################################
############ Directories ############
#####################################
##set directory to all butternut files 
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

########################################
########### Load in Files ##############
########################################
setwd(QUAC_data_files)
##conversion code if needed
QUAC_arp_list <- list.files(path = "QUAC_genind\\garden_wild", pattern = ".arp$")
##now convert 
for(a in 1:length(QUAC_arp_list)){
  
  arp2gen(paste0("QUAC_genind\\garden_wild\\",QUAC_arp_list[[a]]))
  
}

####
QUAC_wild_est_gen <- read.genepop("QUAC_genind\\garden_wild\\QUAC_wild_est.gen", ncode = 3)
QUAC_wild_est_df <- read.csv("QUAC_data_frames\\garden_wild\\QUAC_wild_est_df.csv")
##
QUAC_wild_non_est_gen <- read.genepop("QUAC_genind\\garden_wild\\QUAC_wild_cleaned_non_EST.gen", ncode = 3)
##load in non-est 

##rownames 
rownames(QUAC_wild_est_gen@tab) <- QUAC_wild_est_df$Ind
rownames(QUAC_wild_non_est_gen@tab) <- QUAC_wild_est_df$Ind

##name populations 
levels(QUAC_wild_est_gen@pop) <- unique(QUAC_wild_est_df$Pop)
levels(QUAC_wild_non_est_gen@pop) <- unique(QUAC_wild_est_df$Pop)

########Wild analyses
#########################################
############ Run Analyses ###############
#########################################
##summary
QUAC_est_sum <- summary(QUAC_wild_est_gen)
##run allelic richness 
QUAC_wild_est_alleles <- QUAC_est_sum$pop.n.all/length(QUAC_wild_est_gen@loc.n.all)
QUAC_wild_est_allrich <- colMeans(allelic.richness(QUAC_wild_est_gen)$Ar)	
QUAC_wild_est_allrich_table <- allelic.richness(QUAC_wild_est_gen)$Ar
##now run stats? 
write.csv(QUAC_wild_est_allrich_table,"G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_wild_est_allrich_table.csv")

##read in table 
QUAC_est_org <- read.csv("G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_org.csv")
##no significant difference between populations 
kruskal.test(QUAC_est_org$All_Rich~as.factor(QUAC_est_org$Pop))
##there is a significant difference between loci 
kruskal.test(QUAC_est_org$All_Rich~QUAC_est_org$Locus)
##no difference by pop
boxplot(allelic.richness(QUAC_wild_est_gen)$Ar)
##look at by loci 
boxplot(QUAC_est_org$All_Rich~QUAC_est_org$Locus, names = names(QUAC_wild_est_gen@loc.n.all))

##now try with expected heterozygosity 
QUAC_poppr_EST <- poppr(QUAC_wild_est_gen)

write.csv(QUAC_hexp_df, "G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_est_hexp_df.csv")

##write in hexp 
QUAC_hexp_org <- read.csv("G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_est_hexp_org_df.csv")
##hexp p-value test-- no difference by pop
kruskal.test(QUAC_hexp_org$Hexp~as.factor(QUAC_hexp_org$Pop))
boxplot(QUAC_hexp_org$Hexp~QUAC_hexp_org$Pop)
##hexp 
kruskal.test(QUAC_hexp_org$Hexp~QUAC_hexp_org$Locus)
boxplot(QUAC_hexp_org$Hexp~QUAC_hexp_org$Locus)

#############################################
####### Compare non-est with est ############
#############################################

QUAC_nonest_seppop <- cbind(summary(seppop(QUAC_wild_non_est_gen)$pop1)$Hexp,
                            summary(seppop(QUAC_wild_non_est_gen)$pop2)$Hexp,
                            summary(seppop(QUAC_wild_non_est_gen)$pop3)$Hexp,
                            summary(seppop(QUAC_wild_non_est_gen)$pop4)$Hexp,
                            summary(seppop(QUAC_wild_non_est_gen)$pop5)$Hexp)

write.csv(QUAC_nonest_seppop, "G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_nonest_hexp.csv")
QUAC_hexp_nonest_org <- read.csv("G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_nonest_hexp_org.csv")

##now do a paired t-test
wilcox.test(QUAC_hexp_org$Hexp, QUAC_hexp_nonest_org$Hexp, paired = FALSE, alternative = "two.sided")

##calculate allelic richness 
QUAC_wild_nonest_allrich_table <- allelic.richness(QUAC_wild_non_est_gen)$Ar

##write this out 
write.csv(QUAC_wild_nonest_allrich_table, "G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_wild_non_est_allrich.csv")
QUAC_wild_nonest_allrich <- read.csv("G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_wild_non_est_allrich_org.csv")

###compare between non est and est 
wilcox.test(QUAC_est_org$All_Rich, QUAC_wild_nonest_allrich$All_Rich, alterntive = "one.sided")
pdf("G:\\Shared drives\\Emily_Schumacher\\QUAC_analyses\\QUAC_est_nonest.pdf")
boxplot(QUAC_est_org$All_Rich, QUAC_wild_nonest_allrich$All_Rich, names = c("ESTs", "Non-ESTs"), 
        ylim = c(0,15), col = c("darkseagreen1","darkseagreen4"))
dev.off()
