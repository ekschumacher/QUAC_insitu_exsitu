<b><p><h1 style="color:red;font-size:20px;"> Project Description</b></p></h1>

This repo contains code for determining how well wild genetic diversity for the rare oak species <i>Quercus acerifolia </i> was captured in <i>ex situ</i> collections. (In this repo, <i>Quercus acerifolia</i> is abbreviated to QUAC). <i>Q. acerifolia </i> is a rare oak native to Arkansas with only five wild populations and around 1,000 wild individuals. The wild individuals for this project were collected from 174 wild individuals during a sampling trip by Dr. Sean Hoban and Bailie Fredlock in 2019. The botanic garden individuals were collected in collaboration with 15 different botanic gardens (Bartlett Tree Research Laboratories, Missouri Botanical Garden, The Morton Arboretum, Huntington Botanical Gardens, United States National Arboretum, Arboretum Pouyouleix, Denver Botanic Gardens, Arnold Arboretum, Trompenburg Tuinen & Arboretum, Meise Botanic Garden, Peckerwood, Zoo and BG Plzeň, Morris Arboretum, University of Washington Botanic Gardens, Forstbotanischer Garten Tharandt, Moore Farms Botanical Garden, Chicago Botanic Gardens) and resulted in data collection from 277 botanic garden individuals. 

<p>We aimed to determine if <i>Quercus acerifolia </i> was well conserved in botanic gardens by:</p>
<ol start="1">
<li>Determining if we see a significant difference in the genetic diversity between wild and garden populations</li>
<li>Identifying if all genetic clusters are conserved in botanic gardens</li>
<li>And also identifying what strategies allowed for the best capture of diversity</li>
</ol> 
  
Therefore, to accomplish these goals, we used the code in this repo to clean data genetic, calculate diversity and structure, and then visualize it in figures. 

<b><p><h1 style="color:red;font-size:20px;">Folder Descriptions</b></p></h1>

The two main folders of this project are QUAC_analyses and QUAC_datafiles. 

QUAC_datafiles contains all of the data files used in the analyses of this rare oak. They are broken into several types of data files; QUAC_genind files and QUAC_data_frames. Every genind file has corresponding data frame that stores population and individual names, as well as scores as exported from Geneious. Genind files were all created from Geneious score documents and formatted in Genelax to be loaded in as .arp files which were then converted to adegenet format genind files. Geneclass files were created from the same Genalex files and converted to genepop to be loaded into Geneclass 2. 

QUAC_analyses contains the scripts and results from the analyses. The R Scripts folder contains all of the scripts used for running genetic diversity analyses. The results folder contains all of the raw result files from analyses. 
