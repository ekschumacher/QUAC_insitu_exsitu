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

##load in garden vs. all wild pops
QUAC_garden_allpop_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\QUAC_garden_allwildpop_clean.gen"), ncode = 3)

###rename individuals 
rownames(QUAC_wild_gen@tab) <- QUAC_wild_df$Ind
rownames(QUAC_garden_gen@tab) <- QUAC_garden_df$Ind 

##rename populations 
levels(QUAC_wild_gen@pop) <- c("Porter", "Magazine", "Pryor", "Sugar Loaf", "Kessler")
levels(QUAC_garden_gen@pop) <- unique(QUAC_garden_df$Pop)
levels(QUAC_garden_allpop_gen@pop) <- c("Garden", "Porter", "Magazine", "Pryor", "Sugar Loaf", "Kessler")

###########################
######### PCoA ############
###########################
##now run PCA code 
QUAC_tab <- tab(QUAC_wild_gen, freq=TRUE, NA.method="mean")
##now create PCA tab
QUAC_pca <- dudi.pca(QUAC_tab, center=TRUE, scale=FALSE)
##now look at results
QUAC_pca_df <- cbind(QUAC_pca$li$Axis1, QUAC_pca$li$Axis2)
QUAC_pca_df <- data.frame(QUAC_pca_df)
##add pops and colors 
QUAC_pca_df$Pop <- QUAC_wild_df$Pop
##now add colors 
QUAC_pca_df$colors <- NA
QUAC_pca_df[QUAC_pca_df$Pop == "PorterMountain_W",][,4] <- "cadetblue1"
QUAC_pca_df[QUAC_pca_df$Pop == "MagazineMountain_W",][,4] <- "cadetblue4"
QUAC_pca_df[QUAC_pca_df$Pop == "PryorMountain_W",][,4] <- "darkgreen"
QUAC_pca_df[QUAC_pca_df$Pop == "Sugar_Loaf_Mountains_W",][,4] <- "darkorchid"
QUAC_pca_df[QUAC_pca_df$Pop == "Kessler_Mountain_W",][,4] <- "lightcoral"

##create color vector 
#QUAC_colors <- c("cadetblue1","cadetblue4","darkgreen","darkorchid","lightcoral")
##write out PCoA 
#pdf(paste0(QUAC_analysis_results, "\\QUAC_wildpops_PCA.pdf"))
#s.class(QUAC_pca$li, fac=pop(QUAC_wild_gen),
 #       col=QUAC_colors,
#        axesel=FALSE, cstar=0, cpoint=2, clabel = 1, addaxes = TRUE)
#dev.off()

##percent explained by each axis
pc1 <- signif(((QUAC_pca$eig[1])/sum(QUAC_pca$eig))*100, 3)
pc2 <- signif(((QUAC_pca$eig[2])/sum(QUAC_pca$eig))*100, 3)

##plot PCA
pdf(paste0(QUAC_analysis_results, "\\QUAC_wildpops_better_PCA.pdf"))
plot(QUAC_pca_df$X1, QUAC_pca_df$X2, col = QUAC_pca_df$colors, pch = 17, xlab = paste0("PC1 (", pc1, "%)"), ylab = paste0("PC2 (", pc2, "%)"))
legend('bottomleft', pch = 17, col = c("cadetblue1","cadetblue4","darkgreen","darkorchid","lightcoral"), 
       legend = c("Porter", "Magazine", "Pryor", "Sugar Loaf", "Kessler"))
abline(h = 0)
abline(v = 0)
dev.off()

##now try PCA on garden information 
QUAC_garden_tab <- tab(QUAC_garden_gen, freq=TRUE, NA.method="mean")
##
QUAC_garden_pca <- dudi.pca(QUAC_garden_tab, center=TRUE, scale=FALSE)

##try plotting 
pdf(paste0(QUAC_analysis_results, "\\QUAC_gardenpops_PCA.pdf"))
s.class(QUAC_garden_pca$li, fac=pop(QUAC_garden_gen),
      col= transp(funky(15),.6),
       axesel=FALSE, cstar=0, cpoint=2, clabel = 1, addaxes = TRUE)
dev.off()
##maybe not as useful, try joining everything together and do garden vs. all wild populations 
QUAC_allwildpop_garden_tab <- tab(QUAC_garden_allpop_gen, freq=TRUE, NA.method="mean")
##now run analyses 
QUAC_allwildpop_garden_pca <- dudi.pca(QUAC_allwildpop_garden_tab, center=TRUE, scale=FALSE)

pdf(paste0(QUAC_analysis_results, "\\QUAC_allwildpop_garden_PCA.pdf"))
s.class(QUAC_allwildpop_garden_pca$li, fac=pop(QUAC_garden_allpop_gen),
         col= c("indianred1", "cadetblue1","cadetblue4","darkgreen","darkorchid","lightgoldenrod1"),
        axesel=FALSE, cstar=0, cpoint=2, clabel = 1, addaxes = TRUE)
dev.off()

##calculate percent variation
QUAC_allpopwild_garden_pc1 <- signif(((QUAC_allwildpop_garden_pca$eig[1])/sum(QUAC_allwildpop_garden_pca$eig))*100, 3)
QUAC_allpopwild_garden_pc2 <- signif(((QUAC_allwildpop_garden_pca$eig[2])/sum(QUAC_allwildpop_garden_pca$eig))*100, 3)

