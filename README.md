<b><p><h1 style="color:red;font-size:20px;"> Project Description</b></p></h1>

This repo contains code for determining how well wild genetic diversity for the rare oak species <i>Quercus acerifolia </i> was captured in <i>ex situ</i> collections. (In this repo, <i>Quercus acerifolia</i> is abbreviated to QUAC). <i>Q. acerifolia </i> is a rare oak native to Arkansas with only five wild populations and around 600 wild individuals. Leaf samples of the wild individuals for this project were collected from 174 wild individuals during a sampling trip by colleagues in 2019. Samples of the botanic garden individuals were collected in collaboration with 15 different botanic gardens (Bartlett Tree Research Laboratories, Missouri Botanical Garden, The Morton Arboretum, Huntington Botanical Gardens, United States National Arboretum, Arboretum Pouyouleix, Denver Botanic Gardens, Arnold Arboretum, Trompenburg Tuinen & Arboretum, Meise Botanic Garden, Peckerwood, Zoo and BG Plzeň, Morris Arboretum, University of Washington Botanic Gardens, Forstbotanischer Garten Tharandt, Moore Farms Botanical Garden, Chicago Botanic Gardens) and resulted in data collection from 277 botanic garden individuals. 

<p>We aimed to determine if <i>Quercus acerifolia </i> was well conserved in botanic gardens by:</p>
<ol start="1">
<li>Quantifying the difference in the genetic diversity (using several metrics) between wild and garden populations</li>
<li>Quantifying the percentage of wild occurring alleles which are conserved in botanic gardens</li>
<li>Identifying which wild genetic clusters are conserved in botanic gardens</li>
<li>Identifying what strategies (e.g. minimum sample size) allowed for the best capture of diversity</li>
</ol> 
  
Code in this repo will clean data genetic, calculate diversity and structure, calculate levels of allele capture and create visualizations. 

<b><p><h1 style="color:red;font-size:20px;">Folder Descriptions</b></p></h1>

The two main folders of this project are QUAC_analyses and QUAC_datafiles. 

<b>QUAC_datafiles</b>: This folder contains all of the data files used in the analyses of this rare oak. They are broken into several types of data files; QUAC_genind, QUAC_data_frames, and QUAC_geneclass files. Every genind file has corresponding data frame that stores population and individual names, as well as scores as exported from Geneious. Genind files were all created from Geneious score documents and formatted in Genelax to be loaded in as .arp files which were then converted to adegenet format genind files. Geneclass files were created from the same Genalex files and converted to genepop to be loaded into Geneclass 2 for assignment testing. 

<b>QUAC_analyses</b>: This folder contains all the RScripts used to run analyses on QUAC data files. The Results folder is separated into Clustering analyses (like PCA and STRUCTURE) and Sum_Stats contains all of the overall statistics calculated by wild populations and garden separately. Finally, the Wild_Garden_Comparison folder contains all of the results generated when comparing diversity differences between garden and wild populations. 
