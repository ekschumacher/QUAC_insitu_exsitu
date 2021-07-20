#This code is used to generate PCoAs for the analysis
#This file uses the reduced individual document from the 
#relatedness removal document 

##########################
######## Libraries #######
##########################

library(adegenet)
library(diveRsity)

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
QUAC_garden_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_garden.gen"), ncode = 3)
QUAC_wild_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\QUAC_wild_df.csv"))
QUAC_garden_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\QUAC_garden_df.csv"))

###rename individuals 
rownames(QUAC_wild_gen@tab) <- QUAC_wild_df$Ind
rownames(QUAC_garden_gen@tab) <- QUAC_garden_df$Ind 

##rename populations 
levels(QUAC_wild_gen@pop) <- c("Porter", "Magazine", "Pryor", "Sugar Loaf", "Kessler")
levels(QUAC_garden_gen@pop) <- unique(QUAC_garden_df$Pop)

###########################
######### PCoA ############
###########################

##now run PCA code 
QUAC_tab <- tab(QUAC_wild_gen, freq=TRUE, NA.method="mean")
##now create PCA tab
QUAC_pca <- dudi.pca(QUAC_tab, center=TRUE, scale=FALSE)
##now look at results
s.label(QUAC_pca$li)

s.class(QUAC_pca$li, fac=pop(QUAC_wild_gen), col=funky(15))

##create color vector 
QUAC_colors <- c("cadetblue1","cadetblue4","darkgreen","darkorchid","lightcoral")
##write out 
s.class(QUAC_pca$li, fac=pop(QUAC_wild_gen),
        col=QUAC_colors,
        axesel=FALSE, cstar=0, cpoint=2, clabel = 1)
