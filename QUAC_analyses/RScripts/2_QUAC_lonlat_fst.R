##########################
######## Libraries #######
##########################

library(adegenet)
library(hierfstat)
library(geosphere)

#####################################
############ Directories ############
#####################################
##set directory to all butternut files 
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_diversity\\QUAC_analyses\\Results"

################################
########## Load files ##########
################################
setwd(QUAC_data_files)

QUAC_wild_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_wild.gen"), ncode = 3)
QUAC_wild_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\garden_wild\\QUAC_wild_df.csv"))

###create population name list
QUAC_pop_list <- unique(QUAC_wild_df$Pop)

###rename individuals 
rownames(QUAC_wild_gen@tab) <- QUAC_wild_df$Ind
##name populations
levels(QUAC_wild_gen@pop) <- c("Porter", "Magazine", "Pryor", "Sugar Loaf", "Kessler")

##################################
####### Geographic Analyses ######
##################################
##Calculate mean longitude and latitude for each population
#first create matrices
QUAC_mean_lon <- matrix()
QUAC_mean_lat <- matrix()

##identifying mean latitudes for each population

for(pop in QUAC_pop_list){
  
  QUAC_mean_lon[pop] <- mean(QUAC_wild_df[QUAC_wild_df$Pop == pop,][,3])
  
  
}

for(pop in QUAC_pop_list){
  
  QUAC_mean_lat[pop] <- mean(QUAC_wild_df[QUAC_wild_df$Pop == pop,][,4])
  
  
}

##convert to matrix
QUAC_mean_lon <- matrix(QUAC_mean_lon)
QUAC_mean_lat <- matrix(QUAC_mean_lat)

##document cleanup
QUAC_mean_lon <- QUAC_mean_lon[-1]
QUAC_mean_lat <- QUAC_mean_lat[-1]

#combine into one document for mean long and lat for each pop
QUAC_coords <- matrix(ncol = 2, nrow = length(QUAC_pop_list))
QUAC_coords[,1] <- QUAC_mean_lon
QUAC_coords[,2] <- QUAC_mean_lat
rownames(QUAC_coords) <- c("Porter", "Magazine", "Pryor", "Sugar Loaf", "Kessler")
colnames(QUAC_coords) <- c("Mean Lon", "Mean Lat")

#################################
###### Fst Calculations #########
#################################

##convert to hierfstat 
QUAC_hierfstat <- genind2hierfstat(QUAC_wild_gen)

##run pairwise fst code 
QUAC_fst_df <- pairwise.neifst(QUAC_hierfstat)

##calculate geographic distances between mean locations
QUAC_dist <- matrix(nrow = length(QUAC_pop_list), ncol = length(QUAC_pop_list))

for(first in 1:length(QUAC_pop_list)){
  for(second in 1:length(QUAC_pop_list)){
    QUAC_dist[first,second] <-  distm(QUAC_coords[first,], QUAC_coords[second,], fun = distGeo)/1000
  }
}

##replacce NAs with zeroes 
QUAC_dist[is.na(QUAC_dist)] <- 0

QUAC_fst_df[is.na(QUAC_fst_df)] <- 0

##run relationship
QUAC_fst_dist <- lm(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)])

##plot distance and fst 
pdf(paste0(QUAC_analysis_results, "\\Structure\\QUAC_Dist_Fst.pdf"))
plot(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)], pch = 17, ylim = c(0,0.13), 
     xlim = c(0,200),
     xlab = c("Distance (km)"), ylab = c("Fst"))
abline(QUAC_fst_dist)
legend('bottomleft', legend = c("R2 = -0.12","p-value = 0.865"), bty = 'n')
dev.off()

###test IBD 
#convert to genpop
QUAC_wild_genpop <- genind2genpop(QUAC_wild_gen)
QUAC_gen_dist <- dist.genpop(QUAC_wild_genpop)
QUAC_geo_dist <- dist(QUAC_coords)
ibd <- mantel.randtest(QUAC_gen_dist, QUAC_geo_dist) ##not significant
QUAC_gen_dist_lm <- lm(as.numeric(QUAC_gen_dist)~as.numeric(QUAC_geo_dist))
QUAC_gen_dist_lm_sum <- summary(QUAC_gen_dist_lm)
                       
####
plot(QUAC_gen_dist~QUAC_geo_dist)
abline(QUAC_gen_dist_lm, col="red",lty=2)


