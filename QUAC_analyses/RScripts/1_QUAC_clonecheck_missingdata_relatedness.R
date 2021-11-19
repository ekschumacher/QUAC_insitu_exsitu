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
QUAC_data_files <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_insitu_exsitu\\QUAC_data_files"

QUAC_analysis_results <- "C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_insitu_exsitu\\QUAC_analyses\\Results"

##############################################
############# Convert to Genind ##############
##############################################
##set working directory to load in data files 
setwd(QUAC_data_files)

##convert QUAC all pop document into genind if needed
#arp2gen(paste0(QUAC_data_files, "\\QUAC_genind\\quac_allpop.arp"))

##now read in genind file 
QUAC_allpop_gen <- read.genepop("QUAC_genind\\quac_allpop.gen", ncode = 3)
 
##load files to label population names and indvidual names 
QUAC_allpop_df <- read.csv("QUAC_data_frames\\QUAC_allpop.csv")

##rename individuals in the genind file
rownames(QUAC_allpop_gen@tab) <- QUAC_allpop_df$ID

##rename population names 
QUAC_popnames <- unique(QUAC_allpop_df$POP)
levels(QUAC_allpop_gen@pop) <- QUAC_popnames

##load document for relatedness analysis 
QUAC_rel_df <- read.csv("QUAC_data_frames\\QUAC_rel_df.csv")

#######################################################################
############# Clone check, relatedness, and missing data ##############
#######################################################################
####run clone check 
##convert to clone check format
QUAC_geneclone <- as.genclone(QUAC_allpop_gen)
##identify multi-locus genotypes (non-clones)
QUAC_gen_id <- mlg.id(QUAC_geneclone)
##function to pull out all clone lists
QUAC_clone_index <- which(sapply(QUAC_gen_id,function(x) length(x)>1))
##now remove clones from the matrix 
QUAC_noclones <- clonecorrect(QUAC_geneclone)
##convert back to genind 
QUAC_genind_nocl <- genclone2genind(QUAC_noclones) ##left with 456 individuals 

##remove individuals with too much missing data 
QUAC_genind_nomd <- missingno(QUAC_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##output to genalex file if needed 
#genind2genalex(QUAC_genind_nomd,file="QUAC_data_frames\\QUAC_clone_free_md_free.csv")
##open genalex and conver to arp file to load new genind file 
#back into R, and then also save as a data frame for 
#relatedness analyses

####Relatedness analysis
##run relatedness analysis
QUAC_relatedness_df <- Demerelate(QUAC_rel_df, object = T, value = "loiselle")
##now identify how many individuals have greater than 25% relatedness = half siblings
QUAC_halfsib_names <- names(which(unlist(QUAC_relatedness_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
QUAC_halfsib_names_cleanfront <- gsub("^.*\\.","", QUAC_halfsib_names)

QUAC_halfsib_names_cleanback <- gsub("^.*\\_","", QUAC_halfsib_names_cleanfront)

QUAC_relate_ind_remove <- unique(QUAC_halfsib_names_cleanback) ##260 of these individuals have > 25% relatedness

##now limit genind by these names 
QUAC_relate_red_gen <- QUAC_genind_nomd[!rownames(QUAC_genind_nomd@tab) %in% QUAC_relate_ind_remove,]
##and export to genalex
genind2genalex(QUAC_relate_red_gen, file="QUAC_data_frames\\Relate_Red\\QUAC_relate_red_gen.csv")
##limit data frames 
QUAC_red_relate_df <- QUAC_rel_df[!QUAC_rel_df$Ind %in% QUAC_relate_ind_remove,]
##write out file 
write.csv(QUAC_red_relate_df, "QUAC_data_frames\\Relate_Red\\QUAC_red_relate_df.csv")

