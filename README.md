# QUAC_diversity

Code determining genetic diversity capture of the rare oak species <i>Quercus acerifolia </i> <i>ex situ</i>. (<i>Quercus acerifolia</i> is abbreviated to QUAC in this project). 

The two main folders of this project are QUAC_analyses and QUAC_datafiles. 

QUAC_datafiles contains all of the data files used in the analyses of this rare oak. They are broken into several types of data files; QUAC_genind files and QUAC_data_frames. Every genind file has corresponding data frame that stores population and individual names, as well as scores as exported from Geneious. Genind files were all created from Geneious score documents and formatted in Genelax to be loaded in as .arp files which were then converted to adegenet format genind files. Geneclass files were created from the same Genalex files and converted to genepop to be loaded into Geneclass 2. 

QUAC_analyses contains the scripts and results from the analyses. The R Scripts folder contains all of the scripts used for running genetic diversity analyses. The results folder contains all of the raw result files from analyses. 
