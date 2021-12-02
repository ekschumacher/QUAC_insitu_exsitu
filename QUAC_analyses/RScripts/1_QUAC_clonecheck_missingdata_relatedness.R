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

#####################################
############ Directories ############
#####################################
QUAC_data_files <- "C:/Users/eschumacher/Documents/GitHub/QUAC_insitu_exsitu/QUAC_data_files"

QUAC_analysis_results <- "C:/Users/eschumacher/Documents/GitHub/QUAC_insitu_exsitu/QUAC_analyses/Results"

##############################################
############# Convert to Genind ##############
##############################################
##set working directory to load in data files 
setwd(QUAC_data_files)

##convert QUAC all pop arlequin file into genepop file if needed
#arp2gen(paste0(QUAC_data_files, "/QUAC_genind/quac_allpop.arp"))

##now read in genepop file as a genind for adegenet 
QUAC_allpop_gen <- read.genepop("QUAC_genind/quac_allpop.gen", ncode = 3)
 
##load data frame - used for labeling population and individual names 
QUAC_allpop_df <- read.csv("QUAC_data_frames/QUAC_allpop.csv")

##rename individuals in the genind object 
rownames(QUAC_allpop_gen@tab) <- QUAC_allpop_df$ID

##create population name list and rename the populations in the genind object 
QUAC_popnames <- unique(QUAC_allpop_df$POP)
levels(QUAC_allpop_gen@pop) <- QUAC_popnames

##load relatedness data frame for relatedness analysis 
QUAC_rel_df <- read.csv("QUAC_data_frames/QUAC_rel_df.csv")

#######################################################################
############# Clone check, relatedness, and missing data ##############
#######################################################################
####run clone check 
##convert to genelcone object
QUAC_geneclone <- as.genclone(QUAC_allpop_gen)
##identify multi-locus genotypes (non-clones)
QUAC_gen_id <- mlg.id(QUAC_geneclone)
##function to pull out all clones into a list
QUAC_clone_index <- which(sapply(QUAC_gen_id,function(x) length(x)>1))
##now remove clones from the matrix 
QUAC_noclones <- clonecorrect(QUAC_geneclone)
##convert back to a genind object
QUAC_genind_nocl <- genclone2genind(QUAC_noclones) ##left with 456 individuals 

##remove individuals with too much missing data 
QUAC_genind_nomd <- missingno(QUAC_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##output to genalex file if needed 
#genind2genalex(QUAC_genind_nomd,file="QUAC_data_frames/QUAC_clone_free_md_free.csv")
##open genalex and convert to arp file to load new genepop into R  

###################################
###### Relatedness Analysis #######
###################################
##run relatedness analysis
QUAC_relatedness_df <- Demerelate(QUAC_rel_df, object = T, value = "loiselle")
##now identify how many individuals have greater than 25% relatedness = half siblings
QUAC_halfsib_names <- names(which(unlist(QUAC_relatedness_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
QUAC_halfsib_names_cleanfront <- gsub("^.*\\.","", QUAC_halfsib_names)

QUAC_halfsib_names_cleanback <- gsub("^.*\\_","", QUAC_halfsib_names_cleanfront)

QUAC_relate_ind_remove <- unique(QUAC_halfsib_names_cleanback) ##260 of these individuals have > 25% relatedness

##now limit genind object by these names 
QUAC_relate_red_gen <- QUAC_genind_nomd[!rownames(QUAC_genind_nomd@tab) %in% QUAC_relate_ind_remove,]
##and export to genalex
genind2genalex(QUAC_relate_red_gen, file="QUAC_data_frames/Relate_Red/QUAC_relate_red_gen.csv",
               overwrite = TRUE)
##limit data frame with the indiviudals that are highly related  
QUAC_red_relate_df <- QUAC_rel_df[!QUAC_rel_df$Sample.ID %in% QUAC_relate_ind_remove,]
##write out data frame  
write.csv(QUAC_red_relate_df, "QUAC_data_frames/Relate_Red/QUAC_red_relate_df.csv")

