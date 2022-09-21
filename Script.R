###################################################

###################################################

# Setting work directory

setwd("Documents/Phyllostomidae_Anc_Reconstruction")

# Installing and loading required packages

if(!require("pacman")){
  install.packages("pacman")
}

library(pacman)

pacman::p_load(ape, castor, dplyr, phangorn, phytools)

# Pruned tree and codification matrices were taken from: #
# https://github.com/DanielaP10/Sensibilidad-reconstruccion-de-estados-ancestrales #

# Pruned tree is based on Upham et al. (2019).
