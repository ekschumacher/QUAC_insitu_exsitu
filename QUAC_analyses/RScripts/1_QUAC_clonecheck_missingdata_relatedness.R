###This script details the analyses run on Quercus acerifolia (referred to as QUAC for brevity here) 
#genotype files to prepare them for genetic diversity and structure analyses. 

#########################
#        Libraries      #
#########################

library(diveRsity)
library(adegenet)
library(poppr)
library(Demerelate)

#########################
#   Load Data Files     #
#########################

#set working directory to load in data files 
setwd("../../QUAC_data_files")

##convert QUAC all pop arlequin file into genepop file if needed
#arp2gen("QUAC_adegenet_files/QUAC_allpop.arp")

##now read in genepop file as a genind for adegenet 
QUAC_gen <- read.genepop("QUAC_adegenet_files/QUAC_allpop.gen", ncode = 3)

##load relatedness data frame for relatedness analysis 
QUAC_df <- read.csv("QUAC_data_frames/QUAC_allpop.csv")

##rename individuals in the genind object 
rownames(QUAC_gen@tab) <- QUAC_df$Ind

##create population name list and rename the populations in the genind object 
QUAC_popnames <- unique(QUAC_df$Pop)
levels(QUAC_gen@pop) <- QUAC_popnames

##load in wild data frame of lon/lat coordinates of all individuals 
QUAC_wild_lonlat_df <- read.csv("QUAC_data_frames/QUAC_wild_lonlat.csv")

############################################################
#    Remove Clones and Individuals with missing data       #
############################################################
###clone check 
#convert to genelcone object
QUAC_geneclone <- as.genclone(QUAC_gen)
#identify multi-locus genotypes (non-clones)
QUAC_gen_id <- mlg.id(QUAC_geneclone)
#function to pull out all clones into a list
QUAC_clone_index <- which(sapply(QUAC_gen_id, function(x) length(x)>1))
#create a list of ID with clones
clones_ID <- list()
for(clones in 1:length(QUAC_clone_index)) clones_ID[[clones]] <- QUAC_gen_id[[QUAC_clone_index[[clones]]]]
#now remove clones from the matrix 
QUAC_noclones <- clonecorrect(QUAC_geneclone)
#convert back to a genind object
QUAC_genind_nocl <- genclone2genind(QUAC_noclones) ##left with 456 individuals, there were 7 clones 

##remove individuals with too much missing data 
QUAC_genind_nomd <- missingno(QUAC_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) ##there were 7 individuals with too much missing data, left with 449 individuals

#output to genalex csv if needed 
genind2genalex(QUAC_genind_nomd, file="QUAC_data_frames/Garden_Wild/QUAC_garden_wild_clean_genalex.csv", 
              overwrite = TRUE)

##remove the individuals from the relatedness df that have been removed for missing data and clones
QUAC_clean_df <- QUAC_df[QUAC_df$Ind %in% rownames(QUAC_genind_nomd@tab),]

##write out data frame 
write.csv(QUAC_clean_df, "QUAC_data_frames/Garden_Wild/QUAC_clean_df.csv")

#################################
#      Relatedness Analysis     #
#################################
##Run relatedness analysis on wild and garden individuals separately 
#first separate data frames
QUAC_garden_rel_df <- QUAC_clean_df[QUAC_clean_df$Garden_Wild == "Garden",] ; QUAC_wild_rel_df <- QUAC_clean_df[QUAC_clean_df$Garden_Wild == "Wild",]

##combine files into lists 
QUAC_rel <- list(QUAC_garden_rel_df[,-3], QUAC_wild_rel_df[,-3])
QUAC_pop <- c("garden","wild")

##create half-sib name list 
QUAC_ind_red_list <- list()

##loop to calculate relatedness for garden and wild individuals and reduce data files for relatedness
for(rel in 1:length(QUAC_rel)){
  
  #calculate relatedness between all individuals for either wild or garden df
  QUAC_deme <- Demerelate(QUAC_rel[[rel]], object = T, value = "loiselle")
  
  #extract just the names of highly related inds (half-sib or more)
  QUAC_halfsib_names <- names(which(unlist(QUAC_deme$Empirical_Relatedness) > 0.25))
  
  ##remove extra characters
  QUAC_halfsib_names_cleanfront <-  gsub("^.*\\.","", QUAC_halfsib_names)
  ##remove extra name
  QUAC_wild_halfsib_names_cleanback <- gsub("^.*\\_","", QUAC_halfsib_names_cleanfront)
  
  ##create list of individuals to remove 
  QUAC_ind_red_list[[rel]] <- unique(QUAC_wild_halfsib_names_cleanback)
  
  ##reduce data frame for related individuals 
  QUAC_red_df <- QUAC_rel[[rel]][!QUAC_rel[[rel]]$ID %in% QUAC_ind_red,]
  
  ##write csv of the reduced data frame 
  write.csv(QUAC_red_df, paste0("QUAC_data_frames/Relate_Red/QUAC_", QUAC_pop[[rel]], "_red_df.csv"))
  
}

###reduce genind files 
QUAC_red_gen <- QUAC_genind_nomd[rownames(QUAC_genind_nomd@tab) %in% unlist(list(QUAC_ind_red_list[[1]], QUAC_ind_red_list[[2]])),]

#output to genalex csv if needed 
genind2genalex(QUAC_red_gen, file="QUAC_data_frames/Relate_Red/QUAC_red_genalex.csv", 
               overwrite = TRUE)

