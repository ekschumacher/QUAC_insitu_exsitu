<b><p><h1 style="color:red;font-size:20px;"> Project Description</b></p></h1>

This repo contains code for determining how well wild genetic diversity for the rare oak species <i>Quercus acerifolia </i> was captured in <i>ex situ</i> collections. (In this repo, <i>Quercus acerifolia</i> is abbreviated to QUAC). <i>Q. acerifolia </i> is a rare oak native to Arkansas with only five wild populations and around 600 wild individuals. Leaf samples of the wild individuals for this project were collected from 174 wild individuals during a sampling trips by colleagues in 2019. Samples of the botanic garden individuals were collected in collaboration with 15 different botanic gardens (Bartlett Tree Research Laboratories, Missouri Botanical Garden, The Morton Arboretum, Huntington Botanical Gardens, United States National Arboretum, Arboretum Pouyouleix, Denver Botanic Gardens, Arnold Arboretum, Trompenburg Tuinen & Arboretum, Meise Botanic Garden, Peckerwood, Zoo and BG Plzeň, Morris Arboretum, University of Washington Botanic Gardens, Forstbotanischer Garten Tharandt, Moore Farms Botanical Garden, Chicago Botanic Gardens) and resulted in data collection from 289 botanic garden individuals. 

<p>We aimed to determine if <i>Quercus acerifolia</i> was well conserved in botanic gardens by:</p>
<ol start="1">
<li>Quantifying the difference in the genetic diversity (using several metrics) between wild and garden populations</li>
<li>Quantifying the percentage of wild occurring alleles which are conserved in botanic gardens</li>
<li>Identifying which wild genetic clusters are conserved in botanic gardens</li>
<li>Identifying what strategies (e.g. minimum sample size) allowed for the best capture of diversity</li>
</ol> 
  
Code in this repo cleans genetic data, calculates diversity and structure, and calculates levels of allele capture and create visualizations. 

<b><p><h1 style="color:red;font-size:20px;">Folder Descriptions</b></p></h1>

The two main folders of this project are QUAC_data_files and QUAC_analyses. QUAC_data_files is the folder containing all of the files used in genetic analyses, which are performed in the QUAC_analyses folder. A "cleaned" data file indicates a data files that has individuals with identical reduced to one represenative individual (undergone clone removal) and also individuals with too much missing data removed (25% or more missing data were removed). Starting with 463 individuals, data files were "cleaned" and left with 449 for most genetic analyses. 

<b>QUAC_data_files</b>: This folder contains all of the data files used in the analyses of this rare oak. They are broken into several types of data files; QUAC_adegenet_files, QUAC_data_frames, QUAC_geneclass, and QUAC_structure files. Adegenet files were all created from Geneious score files and formatted in Genalex to be loaded in as .arp files which were then converted to adegenet format genind objects. Geneclass files were created from the same Genalex files and converted to genepop to be loaded into Geneclass 2 for assignment testing. Data frames are loaded in from csv score files that were exported from Geneious. Structure files were created from Genalex files and converted into structure format for loading into the structure program. 

<b>Files within QUAC_data_files</b> 
<ul><li>QUAC_adegenet_files</li></ul>
<ul><ul><li>Description: All files created for genetic analyses in the R package adegenet. "allpop" Genepop files contain all 463 QUAC individuals separated into all 22 populations. "allpop_clean" Genepop files contain 449 QUAC following individual reduction for missing data and clones. These files were then reorganized to create the documents found in the other folders within this folder.</li></ul></ul>
<ul><ul><li>QUAC_allpop.arp</li></ul></ul>
<ul><ul><li>QUAC_allpop.gen</li></ul></ul>
<ul><ul><li>QUAC_allpop_clean.gen</li></ul></ul>
<ul><ul><li>EST</li></ul></ul>
<ul><ul><ul><li>Description: Folder containing Genepop files with loci separated by whether or not they are an expressed sequence tag microsatellite region (EST) or not (non_EST) and also by population type (garden or wild) </li></ul></ul></ul>
<ul><ul><ul><li> QUAC_garden_cleaned_EST_onepop.gen</li></ul></ul></ul>
<ul><ul><ul><li> QUAC_garden_cleaned_onepop_non_EST.gen</li></ul></ul></ul>
<ul><ul><ul><li> QUAC_wild_cleaned_onepop_non_EST.gen</li></ul></ul></ul>

<ul><ul><li>Garden_Wild</li></ul></ul>
<ul><ul><ul><li>Description: Folder containing Genepop file comprising the "cleaned" QUAC individuals separated into only 2 populations - garden and wild. Used for genetic diversity comparison between garden and wild populations. </li></ul></ul></ul>
<ul><ul><ul><li>QUAC_garden_wild_clean.gen</li></ul></ul></ul>

<ul><ul><li>Relate_Red</li></ul></ul>
<ul><ul><ul><li>Description: Folder containing QUAC individuals reduced by relatedness - individuals that were related by more than 25% were reduced to one representative individuals - and then separated by population type - garden or wild. </li></ul></ul></ul>
<ul><ul><ul><li>QUAC_relate_red_garden_wild.gen</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_relate_red_wild.arp</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_relate_red_wild.gen</li></ul></ul></ul>

<b>QUAC_analyses</b>: This folder contains all the RScripts used to run analyses on QUAC data files. The Results folder is separated into Clustering analyses (like PCA and STRUCTURE) and Sum_Stats contains all of the overall statistics calculated by wild populations and garden separately. Finally, the Wild_Garden_Comparison folder contains all of the results generated when comparing diversity differences between garden and wild populations. 
