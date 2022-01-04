##############This script details the analyses run on Quercus aceriflia genotype files to prepare them 
#############for final genetic analyses (referred to as QUAC for brevity here)
###########For the project looking at genetic diversity capture in garden populations compared to wild populations

#########################
#        Libraries      #
#########################

library(diveRsity)
library(adegenet)
library(poppr)
library(Demerelate)

##########################
#    Convert to Genind   #
##########################

#set working directory to load in data files 
setwd("../../QUAC_data_files")

##convert QUAC all pop arlequin file into genepop file if needed
#arp2gen(paste0(QUAC_data_files, "/QUAC_adegenet_files/QUAC_rawscores_allpop.arp"))

##now read in genepop file as a genind for adegenet 
QUAC_garden_wild_gen <- read.genepop("QUAC_adegenet_files/QUAC_rawscores_garden_wild.gen", ncode = 3)

##load relatedness data frame for relatedness analysis 
QUAC_rel_df <- read.csv("QUAC_data_frames/QUAC_rawscores_garden_wild.csv")

##rename individuals in the genind object 
rownames(QUAC_garden_wild_gen@tab) <- QUAC_rel_df$ID

##create population name list and rename the populations in the genind object 
QUAC_popnames <- unique(QUAC_rel_df$POP)
levels(QUAC_garden_wild_gen@pop) <- QUAC_popnames

##load in genepop file with all populations for all individuals 
QUAC_allpop_clean_gen <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_allpop_clean.gen", ncode = 3)
##load in in data file for all pops, cleaned
QUAC_allpop_clean_df <- read.csv("QUAC_data_frames/QUAC_allpop_clean.csv")
#name individuals 
rownames(QUAC_allpop_clean_gen@tab) <- QUAC_allpop_clean_df$Ind
#populations 
QUAC_allpop_names <- unique(QUAC_allpop_clean_df$Pop)
levels(QUAC_allpop_clean_gen@pop) <- QUAC_allpop_names

###################################################
# Remove Clones and Individuals with missing data #
###################################################

###clone check 
#convert to genelcone object
QUAC_geneclone <- as.genclone(QUAC_garden_wild_gen)
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
QUAC_nomd_nocl_df <- QUAC_rel_df[QUAC_rel_df$ID %in% rownames(QUAC_genind_nomd@tab),]

##write out data frame 
write.csv(QUAC_nomd_nocl_rel_df, "QUAC_data_frames/Garden_Wild/QUAC_garden_wild_clean_df.csv")

#########################
# Relatedness Analysis #
########################

##Run relatedness analysis on wild and garden individuals separately 
#first separate data frames
QUAC_garden_rel_df <- QUAC_nomd_nocl_rel_df[QUAC_nomd_nocl_rel_df$POP == "Garden",] ; QUAC_wild_rel_df <- QUAC_nomd_nocl_rel_df[QUAC_nomd_nocl_rel_df$POP == "Wild",]
#first, create wild only genind
QUAC_garden_gen <- seppop(QUAC_genind_nomd)$Garden
QUAC_wild_gen <- seppop(QUAC_genind_nomd)$Wild

##combine files into lists 
QUAC_rel <- list(QUAC_garden_rel_df, QUAC_wild_rel_df)
QUAC_gen <- list(QUAC_garden_gen, QUAC_wild_gen)
QUAC_pop <- c("garden","wild")

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
  QUAC_ind_red <- unique(QUAC_wild_halfsib_names_cleanback)
  
  ##now limit genind by individuals in the half-sibling list 
  QUAC_red_gen <- QUAC_gen[[rel]][!rownames(QUAC_gen[[rel]]@tab) %in% QUAC_ind_red,]
  
  ##write out to Genalex csv
  genind2genalex(QUAC_red_gen, file = paste0("QUAC_data_frames/Relate_Red/QUAC_", QUAC_pop[[rel]], "_genalex.csv"),
                 overwrite = TRUE)
  
  ##reduce data frame for related individuals 
  QUAC_red_df <- QUAC_rel[[rel]][!QUAC_rel[[rel]]$ID %in% QUAC_ind_red,]
  
  ##write csv of the reduced data frame 
  write.csv(QUAC_red_df, paste0("QUAC_data_frames/Relate_Red/QUAC_", QUAC_pop[[rel]], "_df.csv"))
  
}
