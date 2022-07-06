# damping_time
Code used for paper "Reproductive dispersion and damping time scale with life-history speed".

Below are the files needed for analyses and their brief descriptions.
Please either open the R project ("Reproductive dispersion and damping time scale with life-history speed (data and code).Rproj") or edit the file path while running the scripts.

------------------------------------------------------------
Part 1: Data cleanning, calculation and preliminary results
------------------------------------------------------------
> "[1] data processing and plotting for damping paper.R"
 - This is the R code used for data cleanning, calculation and preliminary regression analyses.
 
> "tulja_lab_master_com_p_adre.2.1.RData"
- This is the combined dataset from COMADRE (DOI:10.1111/1365-2656.12482) and COMPADRE (DOI:10.1111/1365-2745.12334) after data correction using original papers.
- See details for modified species/matrices in "COMADRE Mods.pdf".
- Descriptions of dataset can be found in their website (https://compadre-db.org/Data/Comadre).

> "combined dataset of JMG and MO without duplicates.csv"
- This is the GO dataset complied by Jean-Michel Gaillard (DOI: 10.1016/j.tpb.2012.09.003) and Madan Oli (DOI: 10.1016/j.baae.2004.06.002).
- Descriptions of dataset can be found in "Description of JMG and MO dataset.xlsx".

------------------------------------------------------------
Part 2: Figures and Tables in the paper (PGLS, PPCA, etc.)
------------------------------------------------------------
> "[2] phylogenetic analysis for damping paper.R"
- This is the R code used to produce figures in our paper, including results based on PGLS, PPCA, etc.
- Also note that this part of the analyses can be run directly using only the files below.

> "final_tree_corrected_newversion v2.tre"
- This file contains the matched phylogenetic tree using the master phylognetic tree from Open Tree of Life version 12.3 (DOI: 10.5281/zenodo.3937741).

> "match_list_uncorrected_newversion v2.csv"
- This file contains the species that can be matched with the master phylogenetic tree.

> "full animal and plant data v2.csv"
- This file contains the results (life-history traits, such as generation time, reproductive disperson, damping time, etc.) from Part 1.
