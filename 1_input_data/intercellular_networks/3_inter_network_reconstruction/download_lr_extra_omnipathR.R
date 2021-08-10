# Downloading ligand-receptor extra interactions using OmniPathR
# Author: Alberto Valdeolivas

library(OmnipathR)
library(tidyverse)
library(tibble)
library(dplyr)


interactionFormatTransf <- function(InputDf, InteractionType){
  
  OutputInt <- tibble(from = character(), to = character(), 
                      source = character(), database = character())  
  
  n <- nrow(InputDf)
  sources <- dplyr::pull(InputDf, sources)
  sourceNodes <- dplyr::pull(InputDf, from)
  targetNodes <- dplyr::pull(InputDf, to)
  
  for (i in seq(n)){
    currentSources <- unlist(strsplit(sources[i],";"))
    for (j in seq(length(currentSources))){
      OutputInt <- add_row(OutputInt, 
                           from = sourceNodes[i] , 
                           to = targetNodes[i],  
                           # source = paste(currentSources[j], InteractionType, sep="_"),
                           source = currentSources[j],
                           database = currentSources[j]) 
    }
  }
  
  return(OutputInt)
}

lr_Interactions_Omnipath <- import_ligrecextra_interactions() %>%
  dplyr::select(source_genesymbol,target_genesymbol,sources) %>%
  dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
  dplyr::filter(from != to) %>% 
  dplyr::distinct()

## We import Omnipath Inter cellular annotations
InterCell_Annotations <- import_omnipath_intercell()

## We filter those proteins which are mainly annotated as receptor or ligand
Ligands_Receptors <- InterCell_Annotations %>%
  dplyr::filter(category %in% c("receptor","ligand"))

## There are also some complexes. We are going to deal with them by including
## each of its individual proteins in our list
Ligand_Receptors_class <- character()
Ligand_Receptors_name <- character()
for (i in seq(nrow(Ligands_Receptors))){
  if (Ligands_Receptors$entity_type[i] == "complex"){
    Genescomplex <-unlist(strsplit(gsub("COMPLEX:", "", 
                                        Ligands_Receptors$genesymbol[i]),"_"))
    class <- rep(Ligands_Receptors$category[i],length(Genescomplex))
    Ligand_Receptors_name <- c(Ligand_Receptors_name,Genescomplex)
    Ligand_Receptors_class <- c(Ligand_Receptors_class,class)
    
  } else {
    Ligand_Receptors_name <- 
      c(Ligand_Receptors_name, Ligands_Receptors$genesymbol[i]) 
    Ligand_Receptors_class <- 
      c(Ligand_Receptors_class, Ligands_Receptors$category[i]) 
  }
}

## We create a vector with all the ligands and another with all the receptors.
Ligand_Receptors_df <- data.frame(GeneSymbol = Ligand_Receptors_name, 
                                  Class = Ligand_Receptors_class, stringsAsFactors = FALSE) %>%
  dplyr::distinct()
AllLigands_vec <- 
  dplyr::filter(Ligand_Receptors_df, Class == "ligand") %>%
  dplyr::pull(GeneSymbol)
AllReceptors_vec <- 
  dplyr::filter(Ligand_Receptors_df, Class == "receptor") %>%
  dplyr::pull(GeneSymbol)

## We next get protein-protein interactions from the different datasets availabe
## in Omnipath
AllInteractions <- 
  import_post_translational_interactions(exclude = "ligrecextra") %>% 
  dplyr::select(source_genesymbol, target_genesymbol, sources) %>% 
  dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
  dplyr::filter(from != to) %>% 
  dplyr::distinct() 

## I finally match interactions and annotations.
Matching_Interactions_Annotations <- AllInteractions %>%
  dplyr::filter(from %in% AllLigands_vec) %>%
  dplyr::filter(to %in% AllReceptors_vec) %>%
  dplyr::distinct()


##We now combine these two soruces of ligand receptor interactions and then transform to the network format required by the NicheNet method:

## We access to the Omnipath webservice usign the OmnipathR package and we 
## set the number of sources reporting an interactions as its weight 
lr_Network_Omnipath <- 
  bind_rows(lr_Interactions_Omnipath, Matching_Interactions_Annotations) %>%
  dplyr::distinct() %>%
  interactionFormatTransf(InteractionType="LigrecExtra") %>%
  dplyr::distinct() 

## I have to remove self-interactions
lr_Network_Omnipath <- lr_Network_Omnipath %>% 
  dplyr::filter(from != to) 

## I also have to remove interactions going to ligands. See Methods Nichenet 
## paper
ligands <- unique(dplyr::pull(lr_Network_Omnipath, from))
lr_Network_Omnipath <- lr_Network_Omnipath %>% 
  dplyr::filter(!(to %in% ligands))

## There are in addition some records containing not input gene, we remove them
## since they are giving problems with running the model.
lr_Network_Omnipath <- lr_Network_Omnipath %>% 
  dplyr::filter(from != "") %>% 
  dplyr::filter(to != "")
nrow(lr_Network_Omnipath)

saveRDS(lr_Network_Omnipath, 
        'output.rds')

write_tsv(lr_Network_Omnipath, 'output')