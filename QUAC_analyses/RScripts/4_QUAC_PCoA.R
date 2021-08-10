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

########################################
########### Load in Files ##############
########################################
##first, load in reduced genind with no half siblings
QUAC_red_relate_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\Relate_Red\\QUAC_relate_red.gen"), ncode = 3)
##load in data frame 
QUAC_red_relate_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\Relate_Red\\QUAC_relate_red_df.csv"))
###rename individuals 
rownames(QUAC_red_relate_gen@tab) <- QUAC_red_relate_df$Ind
##create population document 
QUAC_pop_names <- unique(QUAC_red_relate_df$Pop)
##name populations in genind file 
levels(QUAC_red_relate_gen@pop) <- QUAC_pop_names

###load files to compare wild and garden individuals 
QUAC_garden_wild_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\Relate_Red\\QUAC_relate_red_garden_wild.gen"), ncode = 3)
##load file to name pops and inds
QUAC_garden_wild_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\Relate_Red\\QUAC_relate_red_garden_wild_df.csv"))
##name inds and pops
rownames(QUAC_garden_wild_gen@tab) <- QUAC_garden_wild_df$Ind
levels(QUAC_garden_wild_gen@pop) <- unique(QUAC_garden_wild_df$Pop)

######load files for all wild pops compared to garden 
QUAC_garden_allwildpop_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\Relate_Red\\QUAC_relate_red_garden_allwildpop.gen"), ncode = 3)
##load data frame 
QUAC_garden_allwildpop_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\Relate_Red\\QUAC_relate_red_garden_allwildpop_df.csv"))
##rename rows and pops 
rownames(QUAC_garden_allwildpop_gen@tab) <- QUAC_garden_allwildpop_df$Ind
levels(QUAC_garden_allwildpop_gen@pop) <- c("Garden","Porter","Magazine","Pryor","Sugarloaf","Kessler")

#####load files for just wild PCA
##genind
QUAC_wildpop_gen <- read.genepop(paste0(QUAC_data_files, "\\QUAC_genind\\Relate_Red\\QUAC_relate_red_allwildpop.gen"), ncode = 3)
##data frame 
QUAC_wildpop_df <- read.csv(paste0(QUAC_data_files, "\\QUAC_data_frames\\Relate_Red\\QUAC_relate_red_allwildpop_df.csv"))
##rename with rows and individual names 
rownames(QUAC_wildpop_gen@tab) <- QUAC_wildpop_df$Ind
levels(QUAC_wildpop_gen@pop) <- c("Porter","Magazine","Pryor","Sugarloaf","Kessler")

###########################
######### PCoA ############
###########################
####first PCA looking at all the different populations 
##now run PCA code 
QUAC_tab <- tab(QUAC_red_relate_gen, freq=TRUE, NA.method="mean")
##now create PCA tab
QUAC_pca <- dudi.pca(QUAC_tab, scale = FALSE)

##create a data frame for the data frame 
QUAC_allpop_pca_df <- data.frame(cbind(QUAC_pca$li$Axis1, QUAC_pca$li$Axis2, QUAC_red_relate_df$Pop))
colnames(QUAC_allpop_pca_df) <- c("Axis1","Axis2","Pop")

##create new color column 
QUAC_allpop_pca_df$Color <- NA

for(i in 1:length(QUAC_pop_names)){
  
   QUAC_allpop_pca_df[QUAC_allpop_pca_df[,3] == QUAC_pop_names[[i]],][,4] <- funky(length(QUAC_pop_names))[[i]]
  
}
##figure out % variation explained 
allpop_pc1 <- signif(((QUAC_pca$eig[1])/sum(QUAC_pca$eig))*100, 3)
allpop_pc2 <- signif(((QUAC_pca$eig[2])/sum(QUAC_pca$eig))*100, 3)
##create pdf of all pops 
pdf(paste0(QUAC_analysis_results, "\\Clustering\\QUAC_allpop_pca.pdf"), width = 10, height = 8)
##create a diagram of all populations 
plot(QUAC_allpop_pca_df[,1], QUAC_allpop_pca_df[,2], pch = 17, col = QUAC_allpop_pca_df[,4], ylim = c(-2,2), xlim = c(-2,3),
     xlab = paste0('PC1 (',allpop_pc1, "%)"), ylab = paste0("PC2 (", allpop_pc2, "%)"))
abline(v = 0)
abline(h = 0)
legend("topright", legend = QUAC_pop_names, border = "black", bty = "o",
       pt.cex = 1, cex = 0.8, pch = 17, col = unique(QUAC_allpop_pca_df[,4]), bg = "white")

dev.off()

##PCA comparing overall wild and garden structure 
#create tab document
QUAC_garden_wild_tab <- tab(QUAC_garden_wild_gen, freq=TRUE, NA.method="mean")
#now run PCA
QUAC_garden_wild_pca <- dudi.pca(QUAC_garden_wild_tab, scale = FALSE, scannf = FALSE, nf = 2)
##create pdf 
pdf(paste0(QUAC_analysis_results, "\\Clustering\\QUAC_garden_wild_pca.pdf"), width = 10, height = 8)
s.class(QUAC_garden_wild_pca$li, fac=pop(QUAC_garden_wild_gen),
        col=c("darkseagreen3", "forestgreen"),
        axesel=FALSE, cstar=0, cpoint=2, clabel = 1, addaxes = TRUE)
dev.off()

###### PCA with garden + all wild pops 
QUAC_garden_allwildpop_tab <- tab(QUAC_garden_allwildpop_gen, freq=TRUE, NA.method="mean")
##run PCA
QUAC_garden_allwildpop_pca <- dudi.pca(QUAC_garden_allwildpop_tab, scannf = 2, scale = FALSE)

##create color vector
QUAC_colors <- c("darkslategray3","dodgerblue2","darkgreen","darkorchid","lightcoral", "cadetblue4")
##write out PCoA 
pdf(paste0(QUAC_analysis_results, "\\Clustering\\QUAC_garden_allwildpop_pca.pdf"))
s.class(QUAC_pca$li, fac=pop(QUAC_garden_allwildpop_gen),
       col=QUAC_colors,
        axesel=FALSE, cstar=0, cpoint=2, clabel = 1, addaxes = TRUE)
dev.off()

#####PCA of all wild populations 
##create tab document
QUAC_wildpop_tab <- tab(QUAC_wildpop_gen, freq=TRUE, NA.method="mean")
##run PCA 
QUAC_wildpop_pca <- dudi.pca(QUAC_wildpop_tab, scannf = FALSE, nf = 2, scale = FALSE)
##calculate % variation explained by each axis 
wildpop_pc1 <- signif(((QUAC_wildpop_pca$eig[1])/sum(QUAC_wildpop_pca$eig))*100, 3)
wildpop_pc2 <- signif(((QUAC_wildpop_pca$eig[2])/sum(QUAC_wildpop_pca$eig))*100, 3)

##write out data file 
pdf(paste0(QUAC_analysis_results, "\\Clustering\\QUAC_wildpop_pca.pdf"), width = 10, height = 8)
s.class(QUAC_wildpop_pca$li, fac = pop(QUAC_wildpop_gen),
        col = QUAC_colors[c(2:6)], 
        axesel=FALSE, cstar=0, cpoint=2, clabel = 1, addaxes = TRUE)
dev.off()
