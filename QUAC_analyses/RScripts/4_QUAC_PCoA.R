###########This code is used to generate PCoAs for the analysis
##########This file uses the reduced individual document from the 
#########relatedness removal document 

##########################
######## Libraries #######
##########################

library(adegenet)
library(diveRsity)

#####################################
########## Load in Files ############
#####################################
setwd("../../QUAC_data_files")

##convert from arlequin files to genepop files if needed
#arp2gen("QUAC_adegenet_files/Relate_Red/QUAC_relate_red_garden_wild.arp")
#arp2gen("QUAC_adegenet_files/Relate_Red/QUAC_relate_red_wildpop.arp")
#arp2gen("QUAC_adegenet_files/Relate_Red/QUAC_relate_red_garden_wildpops.arp")

####Just wild population clustering 
###load reduced genepop wild as a genind object
QUAC_wildpop_gen <- read.genepop("QUAC_adegenet_files/Relate_Red/QUAC_relate_red_wildpop.gen", ncode = 3)
##data frame 
QUAC_wildpop_df <- read.csv("QUAC_data_frames/Relate_Red/QUAC_relate_red_wildpop_df.csv")
##rename genind object with individual names and population names 
rownames(QUAC_wildpop_gen@tab) <- QUAC_wildpop_df$Ind
##create a vector of population name 
QUAC_wildpop_names <- unique(QUAC_wildpop_df$Pop) 
levels(QUAC_wildpop_gen@pop) <- QUAC_wildpop_names

####garden/wild clustering
###load files to compare wild and garden individuals 
QUAC_garden_wild_gen <- read.genepop("QUAC_adegenet_files/Relate_Red/QUAC_relate_red_garden_wild.gen", ncode = 3)
##load file to name pops and inds
QUAC_garden_wild_df <- read.csv("QUAC_data_frames/Relate_Red/QUAC_relate_red_garden_wild_df.csv")
##name inds and pops in genind object
rownames(QUAC_garden_wild_gen@tab) <- QUAC_garden_wild_df$Ind
levels(QUAC_garden_wild_gen@pop) <- unique(QUAC_garden_wild_df$Pop)

###garden in one population compared with all five wild populations 
##load in the genepop file as an adegenet genind object
QUAC_garden_allwildpops_gen <- read.genepop("QUAC_adegenet_files/Relate_Red/QUAC_relate_red_garden_wildpops.gen", ncode = 3)
##now load in the data frame 
QUAC_garden_allwildpops_df <- read.csv("QUAC_data_frames/Relate_Red/QUAC_relate_red_garden_wildpops_df.csv")
##now name individuals and populations within the genind object
rownames(QUAC_garden_allwildpops_gen@tab) <- QUAC_garden_allwildpops_df$Ind
levels(QUAC_garden_allwildpops_gen@pop) <- unique(QUAC_garden_allwildpops_df$Pop)

###########################
######### PCoA ############
###########################
####first run a PCA comparing all of the wild populations 
QUAC_wildpop_tab <- tab(QUAC_wildpop_gen, freq=TRUE, NA.method="mean")
##now create PCA tab
QUAC_wildpop_pca <- dudi.pca(QUAC_wildpop_tab, scale = FALSE)

##create PCA data frame 
QUAC_wildpop_pca_df <- data.frame(cbind(as.numeric(QUAC_wildpop_pca$li$Axis1), as.numeric(QUAC_wildpop_pca$li$Axis2), QUAC_wildpop_df$Pop))
colnames(QUAC_wildpop_pca_df) <- c("Axis1","Axis2","Pop")

##figure out % variation explained 
QUAC_wildpop_pc1 <- signif(((QUAC_wildpop_pca$eig[1])/sum(QUAC_wildpop_pca$eig))*100, 3)
QUAC_wildpop_pc2 <- signif(((QUAC_wildpop_pca$eig[2])/sum(QUAC_wildpop_pca$eig))*100, 3)

##create a ggplot 
setwd("../QUAC_analyses/Results/Clustering")
pdf("QUAC_wildpop_PCA.pdf", width = 10, height = 8)
ggplot(QUAC_wildpop_pca_df, aes(x = as.numeric(QUAC_wildpop_pca_df$Axis1), y = as.numeric(QUAC_wildpop_pca_df$Axis2), 
                                col = Pop)) + xlab(paste0("PC1 (", QUAC_wildpop_pc1, "%)")) + 
  ylab(paste0("PC2 (", QUAC_wildpop_pc2, "%)")) +
  geom_point() + scale_x_continuous( limits = c(-2,2)) + 
  scale_x_continuous(limits = c(-2,2)) +
  scale_color_manual(values=c("royalblue4", "darkolivegreen4", "darkseagreen2",
                            "dodgerblue", "darkorchid")) + theme_bw() + stat_ellipse()

dev.off()

##PCA comparing overall wild and garden structure 
#create tab document
QUAC_garden_wild_tab <- tab(QUAC_garden_wild_gen, freq=TRUE, NA.method="mean")
#now run PCA
QUAC_garden_wild_pca <- dudi.pca(QUAC_garden_wild_tab, scale = FALSE, scannf = FALSE, nf = 2)

##calculate % variation explained by axis 
QUAC_garden_wild_pc1 <- signif(((QUAC_garden_wild_pca$eig[1])/sum(QUAC_garden_wild_pca$eig))*100, 3)
QUAC_garden_wild_pc2 <- signif(((QUAC_garden_wild_pca$eig[2])/sum(QUAC_garden_wild_pca$eig))*100, 3)

##create pdf 
pdf("QUAC_garden_wild_pca.pdf", width = 10, height = 8)
s.class(QUAC_garden_wild_pca$li, fac=pop(QUAC_garden_wild_gen),
        col=c("darkseagreen3", "forestgreen"),
        axesel=FALSE, cstar=0, cpoint=2, clabel = 1, addaxes = TRUE)
dev.off()

#####PCA of garden vs. all wild populations populations 
##create tab document
QUAC_garden_allwildpops_tab <- tab(QUAC_garden_allwildpops_gen, freq=TRUE, NA.method="mean")
##run PCA 
QUAC_garden_allwildpops_pca <- dudi.pca(QUAC_garden_allwildpops_tab, scannf = FALSE, nf = 2, scale = FALSE)
##calculate % variation explained by each axis 
QUAC_garden_allwildpops_pc1 <- signif(((QUAC_garden_allwildpops_pca$eig[1])/sum(QUAC_garden_allwildpops_pca$eig))*100, 3)
QUAC_garden_allwildpops_pc2 <- signif(((QUAC_garden_allwildpops_pca$eig[2])/sum(QUAC_garden_allwildpops_pca$eig))*100, 3)

##create a data frame 
QUAC_garden_allwildpops_pca_df <- data.frame(cbind(QUAC_garden_allwildpops_pca$li$Axis1,
                                                   QUAC_garden_allwildpops_pca$li$Axis2,
                                                   QUAC_garden_allwildpops_df$Pop))

##put in colnames
colnames(QUAC_garden_allwildpops_pca_df) <- c("Axis1","Axis2","Pop")

##write out data file 
pdf("QUAC_garden_wildpops_pca.pdf", width = 10, height = 8)

ggplot(QUAC_garden_allwildpops_pca_df, aes(x = as.numeric(QUAC_garden_allwildpops_pca_df$Axis1), 
                                           y = as.numeric(QUAC_garden_allwildpops_pca_df$Axis2), 
                                col = Pop)) + xlab(paste0("PC1 (", QUAC_garden_allwildpops_pc1, "%)")) + 
  ylab(paste0("PC2 (", QUAC_garden_allwildpops_pc2, "%)")) +
  geom_point() + scale_x_continuous( limits = c(-2,2)) + 
  scale_x_continuous(limits = c(-2,2)) +
  scale_color_manual(values=c("hotpink","royalblue4", "darkolivegreen4", "darkseagreen2",
                              "dodgerblue", "darkorchid")) + theme_bw() + stat_ellipse()

dev.off()

sessionInfo()
