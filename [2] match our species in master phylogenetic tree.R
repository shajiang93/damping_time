## =====================================
## 2, match species from master phylogenetic tree ####
## =====================================
## it takes a very long time to run
## the master phylogenetic tree is too large to upload in github
## for phylogenetic analysis results in the paper, check file 3
install.packages("remotes")
remotes::install_github("moshagen/multiTreeR")
library(tidyverse)
library(mulTree) # need to find a way to download these packages
library(treeio)

rm(list=ls())
setwd("/Users/shajiang/Documents/AAA @Stanford/research/Notes_AgeStructure")
getwd()

source("./phylogenetic tree/Demography_functions.R")
source("./phylogenetic tree/phylo_bind_functions.R")

#### use Open life history Tree ####
met_tree<-read.tree(file="./phylogenetic tree/labelled_supertree_ottnames.tre")
Ntip(met_tree)
met_tree  <- makeLabel(met_tree) 
all_data <- read.csv2("./output data/full animal and plant data v2.csv")

##==== get our list of species that we will build the tree for ==== ####
species_list <- data.frame(species =(all_data$SpeciesAccepted), 
                           class = (all_data$Class), 
                           phyla = (all_data$Phylum),
                           stringsAsFactors=FALSE)
species_clean_list_u <- unique(species_list$species)

full_uni <- data.frame(species_grep = species_clean_list_u, 
                       species = species_clean_list_u , 
                       species_match = rep(0, length(unique(species_clean_list_u))), 
                       matched = rep(0, length(unique(species_clean_list_u))), 
                       taxonomic_class = rep(0, length(unique(species_clean_list_u))),
                       species.original = species_clean_list_u,
                       stringsAsFactors=FALSE)

### revise species name in the unmatched data
species_name_correct <- as.vector(full_uni$species_grep)
species_name_correct[species_name_correct == "Taraxacum campylodes"] <- "Taraxacum officinale"
species_name_correct[species_name_correct == "Carya sinensis"] <- "Annamocarya sinensis"
species_name_correct[species_name_correct == "Platynathrus peltifer"] <- "Platynothrus peltifer"
species_name_correct[species_name_correct == "Cornu aspersa"] <- "Helix aspersa"
species_name_correct[species_name_correct == "Histrix africaeaustralis"] <- "Hystrix africaeaustralis"
species_name_correct[species_name_correct == "Lutra canadensis"] <- "Lontra canadensis"
species_name_correct[species_name_correct == "Mirounga Leonina"] <- "Mirounga leonina"
species_name_correct[species_name_correct == "Pathera leo"] <- "Panthera leo"
species_name_correct[species_name_correct == "Zapus hodsonius"] <- "Zapus hudsonius"

full_uni$species_grep <- species_name_correct
full_uni$species <- species_name_correct

###remove all the subspecies information
full_uni[,1]  <- gsub(" subsp.*","", full_uni[,1])
full_uni[,1]  <- gsub(" ","_", full_uni[,1])

species_list[,1]  <- gsub(" subsp.*","", species_list[,1])
species_list[,1]  <- gsub(" ","_", species_list[,1])

# have duplicate data after name correction
full_uni$species_grep[duplicated(full_uni$species_grep)==T]
full_uni <- full_uni[!duplicated(full_uni$species_grep), ] 

#This pulls out all the entries in the open tree of life phylogeny that correspond 
#with the species names in our COMADRE species list.

for(i in 1:( length(full_uni[,1]))){
  if(any(grep(full_uni[i,1], met_tree$tip.label)) == T){
    full_uni[i,3]  <- met_tree$tip.label[grep(full_uni[i,1], met_tree$tip.label)][1]
    full_uni[i,4]  <- "yes"
  }
  
  else{
    full_uni[i,3]  <- full_uni[i,1]
    full_uni[i,4]  <- "no"
    full_uni[i,5] <- species_list[species_list$speciesAcc == c(full_uni[i,2]),"class"][1]
  }
  print(i)
}

###addition of Class information
for(i in 1:( length(full_uni[,1]))){
  
  full_uni[i,5] <- as.vector(species_list[(grep(full_uni[i,1],species_list[,1]))[1],"class"])
}
full_uni%>%
  mutate(name.same = as.vector(full_uni$species %in% full_uni$species.original))
full_uni%>%
  group_by(matched)%>%
  count()

length(unique(full_uni$species_match))

## write tree ####
length(unique(met_tree$tip.label))
length(unique(met_tree$node.label))

final_tree <- comparative.data(phy = met_tree, 
                               data = full_uni, 
                               names.col = "species_match", 
                               force.root = TRUE)$phy
plot(final_tree, cex = 0.2, type = "fan")

for(i in 1:(length(final_tree$tip.label))){
  
  if(any(final_tree$tip.label[i] == full_uni$species_match)){
    
    final_tree$tip.label[i] <-  gsub("_", 
                                     " ",
                                     full_uni[final_tree$tip.label[i] ==
                                                full_uni$species_match,
                                              "species"])
  } 
}


final_tree_length <-compute.brlen(final_tree, 
                                  method = "Grafen",
                                  power = 1)

write.tree(phy = final_tree_length,
           file = "./phylogenetic tree/final_tree_uncorrected_newversion v2.tre")
write.table(full_uni, file="./phylogenetic tree/match_list_uncorrected_newversion v2.csv",
            sep = ",",row.names = T)

