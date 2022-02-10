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

<ul><li>QUAC_data_frames</li></ul>
<ul><ul><li>Description: This folder contains all data frames created for genetic analyses. "cleaned" data frames reflect that individuals have been removed for missing data and clones. CSV refers to a comma separated value file.</li></ul></ul>
<ul><ul><ul><li>QUAC_garden_occ.csv: A CSV of all the garden individuals with source GPS locations</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wild_lon_lat.csv: CSV of QUAC wild individuals and their sample location </li></ul></ul></ul>
<ul><ul><ul><li>QUAC_allpop.csv: CSV of QUAC wild individuals and their sample location </li></ul></ul></ul>

<ul><ul><li>EST</li></ul></ul>
<ul><ul><ul><li>Description: Folder containing QUAC individuals separated by population type - garden or botanic garden - with loci separated by EST vs. genomic microsatellites.</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_garden_clean_EST_df.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_garden_clean_non_EST_df.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wild_clean_EST_df.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wild_clean_non_EST_df.csv</li></ul></ul></ul>

<ul><ul><li>Garden_Wild</li></ul></ul>
<ul><ul><ul><li>Description: Folder containing QUAC data frames created with populations and population type for each individual. </li></ul></ul></ul>
<ul><ul><ul><li>QUAC_clean_df.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_garden_wild_clean_genalex.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wild_lonlat_allpop_clean_df.csv</li></ul></ul></ul>

<ul><ul><li>Relate_Red</li></ul></ul>
<ul><ul><ul><li>Description: Folder containing QUAC data frames created by relatedness reduction - individuals with 25% missing data are removed - that are used for clustering analyses.</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_garden_relate_red_df.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_garden_relate_red_genalex.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wild_relate_red_df.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wild_relate_red_genalex.csv</li></ul></ul></ul>

<ul><li>QUAC_geneclass</li></ul>
<ul><ul><li>Description: This folder contains all data frames created for assignment testing in Geneclass. Input files are in a Genepop format and were generated before and after reduction for relatedness (individuals related by 25% or more were removed to created "relate_red" data files). CSV refers to a comma separated value file.</li></ul></ul>
<ul><ul><ul><li>QUAC_garden_input.genepop</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_garden_relate_red_input.genepop</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_geneclass_output.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_geneclass_relate_red_output.csv</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wild_relate_red_input.genepop</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wildpops_woKessler.genepop</li></ul></ul></ul>


<ul><li>QUAC_strucure</li></ul>
<ul><ul><li>Description: This folder contains text files generated from relatedness reduced Genalex files used to generate structure diagrams and so were exported from Genalex in structure format. Only "READY" files were used in structure analyses, either "str" only document were used as conversion files. </li></ul></ul>
<ul><ul><ul><li>QUAC_str.txt</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_str_READY.txt</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wildonly_str.txt</li></ul></ul></ul>
<ul><ul><ul><li>QUAC_wildonly_str_READY.txt</li></ul></ul></ul>

<b>QUAC_analyses</b>: This folder contains all the RScripts used to run analyses on QUAC data files. The Results folder is separated into Clustering analyses (like PCA and STRUCTURE) and Sum_Stats contains all of the overall statistics calculated by wild populations and garden separately. Finally, the Wild_Garden_Comparison folder contains all of the results generated when comparing diversity differences between garden and wild populations. 

<ul><li>Results</li></ul>
<ul><ul><li>Description: This folder contains the results of the analyses performed in the RScripts folder.</li></ul></ul>

<ul><ul><ul><li>Clustering: This folder contains the results of clustering analyses; broadly including PCA results, structure results, and distance/fst linear relationships, and the results of Geneclass assignment testing.</li></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_Dist_Fst.jpg</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_Dist_Fst.pdf</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_geneclass.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_geneclass_relate_red.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_PCA_garden_allwildpops.pdf</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_PCA_garden_wild.pdf</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_PCA_wild.pdf</li></ul></ul></ul></ul>

<ul><ul><ul><li>Sum_Stats: This folder contains the general results of genetic summary stats analyses - calculating  allelic richness and heterozygosity and multilocus genotype within all botanic gardens and wild populations.</li></ul></ul></ul>
<ul><ul><ul><ul><li>Garden_gendiv_sumstat_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>HWE_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>LD_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>Null_All_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>Wild_gendiv_sumstat_df.csv</li></ul></ul></ul></ul>

<ul><ul><ul><li>Wild_Garden_Comparison: This folder contains the results of comparisons between population types - botanic garden individuals as a whole group and wild individuals as a group. The significance comparison is stored here between genetic diversity metrics and the comparisons between EST and genomic microsatellites, but also percent alleles captured and calculations of minimum sample sizes to capture diversity are stored in here. There are 2 types of files in several cases - ndrop0 and ndrop2. Alleles with 2 or fewer frequencies are sometimes dropped (ndrop2) because they may be genotyping mistakes or very rare alleles. Therefore, we performed analyses both with and without these alleles to determine how results differed with and without these alleles. </li></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_all_cap_garden_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_all_cap_garden_df_dr_0.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_all_exist_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_all_exist_df_dr_0.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_all_resampling_ndrop0.pdf</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_all_resampling_ndrop2.pdf</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_EST_boxplots.pdf</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_EST_mean_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_gendiv_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_min_samp_95.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_wild_cap_df.csv</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC_wild_cap_df_dr_0.csv</li></ul></ul></ul></ul>
