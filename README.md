# damping_time
Code used for paper "Reproductive dispersion and damping time scale with life-history speed"

------------------------------------------------------------
Part 1: Data cleanning, calculation and preliminary results
------------------------------------------------------------
Files needed:
"[1] data processing and plotting for damping paper.R"
 - This is the R code used for data cleanning, calculation and preliminary regression analyses
 
"tulja_lab_master_com_p_adre.2.1.RData"
- This is the combined dataset from COMADRE (DOI:10.1111/1365-2656.12482) and COMPADRE (DOI:10.1111/1365-2745.12334) after data correction using original papers.
- See details for modified species/matrices in "COMADRE Mods.pdf".

"combined dataset of JMG and MO without duplicates.csv"
- This is the GO dataset complied by Jean-Michel Gaillard (DOI: 10.1016/j.tpb.2012.09.003) and Madan Oli (DOI: 10.1016/j.baae.2004.06.002).

------------------------------------------------------------
Part 2: Phylogenetic Tree
-------------------------------------------------------------
"[2] match our species in master phylogenetic tree.R"
- This is the R code used to match our data with the phylognetic tree from Open Tree of Life version 12.3 (DOI: 10.5281/zenodo.3937741)
- The above code requires two functions: "Demography_functions.R" and "phylo_bind_functions.R".

"labelled_supertree_ottnames.tre.zip"
- This is the phylognetic tree from Open Tree of Life version 12.3 (DOI: 10.5281/zenodo.3937741). Please unzip before using.

------------------------------------------------------------
Part 3: Results in the paper (PGLS, PPCA, etc.)
------------------------------------------------------------
"[3] phylogenetic analysis for damping paper.R"
- This is the R code used to produce figures in our paper, including results based on PGLS, PPCA, etc.
