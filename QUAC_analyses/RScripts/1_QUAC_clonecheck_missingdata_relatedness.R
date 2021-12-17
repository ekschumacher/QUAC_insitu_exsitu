##############This script details the preliminary analyses
#############conducted on Quercus acerifolia 
############(referred to as QUAC for brevity here)
###########For the project looking at genetic diversity capture 
##########in garden populations compared to wild populations

##########################
######## Libraries #######
##########################

library(diveRsity)
library(adegenet)
library(poppr)
library(Demerelate)

##############################################
############# Convert to Genind ##############
##############################################
##set working directory to load in data files 
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

###load in all population files 
##load in genepop file with all populations for all individuals 
QUAC_allpop_clean_gen <- read.genepop("QUAC_adegenet_files/Garden_Wild/QUAC_allpop_clean.gen", ncode = 3)
##load in in data file for all pops, cleaned
QUAC_allpop_clean_df <- read.csv("QUAC_data_frames/QUAC_allpop_clean.csv")
#name individuals 
rownames(QUAC_allpop_clean_gen@tab) <- QUAC_allpop_clean_df$Ind
#populations 
QUAC_allpop_names <- unique(QUAC_allpop_clean_df$Pop)
levels(QUAC_allpop_clean_gen@pop) <- QUAC_allpop_names

#######################################################################
############# Clone check, relatedness, and missing data ##############
#######################################################################
####run clone check 
##convert to genelcone object
QUAC_geneclone <- as.genclone(QUAC_garden_wild_gen)
##identify multi-locus genotypes (non-clones)
QUAC_gen_id <- mlg.id(QUAC_geneclone)
##function to pull out all clones into a list
QUAC_clone_index <- which(sapply(QUAC_gen_id, function(x) length(x)>1))
##create a list of ID with clones
clones_ID <- list()
for(clones in 1:length(QUAC_clone_index)) clones_ID[[clones]] <- QUAC_gen_id[[QUAC_clone_index[[clones]]]]
##now remove clones from the matrix 
QUAC_noclones <- clonecorrect(QUAC_geneclone)
##convert back to a genind object
QUAC_genind_nocl <- genclone2genind(QUAC_noclones) ##left with 456 individuals, there were 7 clones 

##remove individuals with too much missing data 
QUAC_genind_nomd <- missingno(QUAC_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) ##there were 7 individuals with too much missing data, left with 449 individuals

##output to genalex file if needed 
genind2genalex(QUAC_genind_nomd, file="QUAC_data_frames/Garden_Wild/QUAC_garden_wild_clean_genalex.csv", 
              overwrite = TRUE)
##open genalex and convert to arp file to load new genepop into R  
##remove the individuals from the relatedness df that have been removed for missing data and clones
QUAC_nomd_nocl_rel_df <- QUAC_rel_df[QUAC_rel_df$ID %in% rownames(QUAC_genind_nomd@tab),]

##write out data frame 
write.csv(QUAC_nomd_nocl_rel_df, "QUAC_data_frames/Garden_Wild/QUAC_garden_wild_clean_df.csv")

###################################
###### Relatedness Analysis #######
###################################
###########Run relatedness analysis on wild and garden individuals separately 
###first separate data frames
QUAC_garden_rel_df <- QUAC_nomd_nocl_rel_df[QUAC_nomd_nocl_rel_df$POP == "Garden",]
QUAC_wild_rel_df <- QUAC_nomd_nocl_rel_df[QUAC_nomd_nocl_rel_df$POP == "Wild",]

##write out data frames
write.csv(QUAC_garden_rel_df, "QUAC_data_frames/Garden_Wild/QUAC_garden_clean_df.csv")
write.csv(QUAC_wild_rel_df, "QUAC_data_frames/Garden_Wild/QUAC_wild_clean_df.csv")

#########wild relatedness
##run relatedness analysis on the cleaned score data frame 
QUAC_wild_relatedness_df <- Demerelate(QUAC_wild_rel_df, object = T, value = "loiselle")
##now identify how many individuals have greater than 25% relatedness = half siblings
QUAC_wild_halfsib_names <- names(which(unlist(QUAC_wild_relatedness_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
QUAC_wild_halfsib_names_cleanfront <- gsub("^.*\\.","", QUAC_wild_halfsib_names)

QUAC_wild_halfsib_names_cleanback <- gsub("^.*\\_","", QUAC_wild_halfsib_names_cleanfront)

QUAC_wild_relate_ind_remove <- unique(QUAC_wild_halfsib_names_cleanback)

##now limit genind 
#first, create wild only genind
QUAC_garden_nomd_nocl_gen <- seppop(QUAC_genind_nomd)$Garden
QUAC_wild_nomd_nocl_gen <- seppop(QUAC_genind_nomd)$Wild
#now limit by the highly related individuals
QUAC_wild_relate_red_gen <- QUAC_wild_nomd_nocl_gen[!rownames(QUAC_wild_nomd_nocl_gen@tab) %in% QUAC_wild_relate_ind_remove,]
##export genind object as a genalex file 
genind2genalex(QUAC_wild_relate_red_gen, file="QUAC_data_frames/Relate_Red/QUAC_relate_red_wild_genalex.csv",
               overwrite = TRUE)
##limit data frame with the indiviudals that are highly related  
QUAC_red_relate_garden_wild_df <- QUAC_rel_df[!QUAC_rel_df$ID %in% QUAC_relate_ind_remove,]
##write out data frame  
write.csv(QUAC_red_relate_garden_wild_df, "QUAC_data_frames/Relate_Red/QUAC_red_relate_garden_wild_df.csv")

#####prep all wild pop file for structure  
##reduce genind object for only wild populations 
QUAC_wild_allpop_clean_gen <- QUAC_allpop_clean_gen[rownames(QUAC_allpop_clean_gen@tab) %in% rownames(QUAC_wild_nomd_nocl_gen@tab),]
##reduce genind 
QUAC_wild_relate_red_allpop_gen <- QUAC_wild_allpop_clean_gen[!rownames(QUAC_wild_allpop_clean_gen@tab) %in% QUAC_wild_relate_ind_remove,]
##export genind object as a genalex file 
genind2genalex(QUAC_wild_relate_red_allpop_gen, file="QUAC_data_frames/Relate_Red/QUAC_wild_relate_red_allpop_genalex.csv",
               overwrite = TRUE)

########limiting garden individuals by relatedness
QUAC_garden_relatedness_df <- Demerelate(QUAC_garden_rel_df, object = T, value = "loiselle")
##now identify how many individuals have greater than 25% relatedness = half siblings
QUAC_garden_halfsib_names <- names(which(unlist(QUAC_garden_relatedness_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
QUAC_garden_halfsib_names_cleanfront <- gsub("^.*\\.","", QUAC_garden_halfsib_names)

QUAC_garden_halfsib_names_cleanback <- gsub("^.*\\_","", QUAC_garden_halfsib_names_cleanfront)

QUAC_garden_relate_ind_remove <- unique(QUAC_garden_halfsib_names_cleanback)

##reduce genind object for only garden populations 
QUAC_garden_allpop_clean_gen <- QUAC_allpop_clean_gen[rownames(QUAC_allpop_clean_gen@tab) %in% rownames(QUAC_garden_nomd_nocl_gen@tab),]
##now reduce by relatedness
QUAC_garden_relate_red_allpop_gen <- QUAC_garden_allpop_clean_gen[!rownames(QUAC_garden_allpop_clean_gen@tab) %in% QUAC_garden_relate_ind_remove,]
##export to genalex file 
genind2genalex(QUAC_garden_relate_red_allpop_gen, file="QUAC_data_frames/Relate_Red/QUAC_garden_relate_red_allpop_genalex.csv",
               overwrite = TRUE)

