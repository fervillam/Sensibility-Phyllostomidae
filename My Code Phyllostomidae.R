
# Installing and loading required packages

if(!require("pacman")){
  install.packages("pacman")
}

library(pacman)

pacman::p_load(ape, castor, dplyr, phytools)

# Building the example topology that will be used in the analysis. This tree
# was taken from Paradis (2012)

Primates.Tree <- read.tree("Documents/Comparada II/Paradis 2012/primfive.tre")
summary(Primates.Tree) # It has five terminals
plot(Primates.Tree)

# Creating different transformation matrices

  # One parameter: equal rates

  Equal.Two <- matrix(c(0, 1, 1, 0), nrow = 2) 
  
  # Two parameters: unequal rates
  
  Unequal.Two <- matrix(c(0, 1, 2, 0), nrow = 2)
  
  # Three states
  
    # Symmetrical and equal rates
  
    Sym.Equal.Three <- matrix(c(0, rep(1, 3), 0, rep(1, 3),  0), nrow = 3)
    
    # Symmetrical and unequal rates
    
    Sym.Unequal.Three <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3)
    
    # Asymmetrical and unequal rates
    
    Asym.Unequal.Three <- matrix(c(0, 1:3, 0, 4:6, 0), nrow = 3)

    # Cyclical
    
    Cycl.Three <- matrix(c(0, 0, 3, 1, 0, 0, 0, 2, 0), nrow = 3)
    
# Ancestral reconstruction of discrete characters: nostrils
    
  # Maximum Likelihood #
    
  # Binary matrices: Equal Rates
    
    # Rhinarium
    
    Primates.Rhinarium <- c(1, 1, 1, 1, 2)
    
    Rhinarium.Anc.Rec.ER.Root <- ace(Primates.Rhinarium, Primates.Tree, type = "d",
                             model = "ER")
    
    Rhinarium.Anc.Rec.ER.Nodes <- Rhinarium.Anc.Rec.ER.Root$lik.anc
    
    # Flat nostrils
    
    Primates.Flat.Nost <- c(1, 1, 1, 2, 1)
    
    Flat.Nost.Anc.Rec.ER.Root <- ace(Primates.Flat.Nost, Primates.Tree, type = "d",
                                  model = "ER")
    
    Flat.Nost.Anc.Rec.ER.Nodes <- Flat.Nost.Anc.Rec.ER.Root$lik.anc
    
    # Downward nostrils
    
    Primates.Down.Nost <- c(2, 2, 2, 1, 1)
    
    Down.Nost.Anc.Rec.ER.Root <- ace(Primates.Down.Nost, Primates.Tree, type = "d",
                                  model = "ER")
    
    Down.Nost.Anc.Rec.ER.Nodes <- Down.Nost.Anc.Rec.ER.Root$lik.anc
    
  # Binary matrices: All-Rates-Different
  
    # Rhinarium
    
    Rhinarium.Anc.Rec.ARD.Root <- ace(Primates.Rhinarium, Primates.Tree, type = "d",
                                      model = Unequal.Two)
    
    Rhinarium.Anc.Rec.ARD.Nodes <- Rhinarium.Anc.Rec.ARD.Root$lik.anc
    
    # Flat nostrils
    
    Flat.Nost.Anc.Rec.ARD.Root <- ace(Primates.Flat.Nost, Primates.Tree, type = "d",
                                      model = Unequal.Two)
    
    Flat.Nost.Anc.Rec.ARD.Nodes <- Flat.Nost.Anc.Rec.ARD.Root$lik.anc
    
    # Downward nostrils
    
    Down.Nost.Anc.Rec.ARD.Root <- ace(Primates.Down.Nost, Primates.Tree, type = "d",
                                      model = Unequal.Two)
    
    Down.Nost.Anc.Rec.ARD.Nodes <- Down.Nost.Anc.Rec.ARD.Root$lik.anc
    
  # Multistate matrix: Equal Rates 
    
    # where rhinarium is "1", flat nostril "2" and downward nostrils "3" #
    
  Primates.Mult.Nost <- c(3, 3, 3, 2, 1)
    
  Mult.Anc.Rec.ER.Root <- ace(Primates.Mult.Nost, Primates.Tree, type = "d",
                                model = "ER")
    
  Mult.Anc.Rec.ER.Nodes <- Mult.Anc.Rec.ER.Root$lik.anc
    
  # Multistate matrix: Unequal Rates
    
  Mult.Anc.Rec.ARD.Root <- ace(Primates.Mult.Nost, Primates.Tree, type = "d",
                               model = Sym.Unequal.Three)
  
  Mult.Anc.Rec.ARD.Nodes <- Mult.Anc.Rec.ARD.Root$lik.anc
  
  
  # Maximum Parsimony #
  
  # Binary matrices: Equal Rates
  
    # Rhinarium
  
    Rhinarium.Pars.Rec.ER <- asr_max_parsimony(Primates.Tree, Primates.Rhinarium,
                                               transition_costs = "all_equal",
                                               Nstates = 2)
    
    # Flat nostrils
    
    Flat.Nost.Pars.Rec.ER <- asr_max_parsimony(Primates.Tree, Primates.Flat.Nost,
                                               transition_costs = "all_equal",
                                               Nstates = 2)
    
    # Downward nostrils
    
    Down.Nost.Pars.Rec.ER <- asr_max_parsimony(Primates.Tree, Primates.Down.Nost,
                                               transition_costs = "all_equal",
                                               Nstates = 2)
    
  # Binary matrices: Unequal Rates
    
    # Rhinarium
    
    Rhinarium.Pars.Rec.UR <- asr_max_parsimony(Primates.Tree, Primates.Rhinarium,
                                               transition_costs = Unequal.Two,
                                               Nstates = 2)
    
    # Flat nostrils
    
    Flat.Nost.Pars.Rec.UR <- asr_max_parsimony(Primates.Tree, Primates.Flat.Nost,
                                               transition_costs = Unequal.Two,
                                               Nstates = 2)
    
    # Downward nostrils
    
    Down.Nost.Pars.Rec.UR <- asr_max_parsimony(Primates.Tree, Primates.Down.Nost,
                                               transition_costs = Unequal.Two,
                                               Nstates = 2)
    
  # Multistate matrix: Equal Rates
    
  Mult.Pars.Rec.ER <- asr_max_parsimony(Primates.Tree, Primates.Mult.Nost,
                                        transition_costs = "all_equal",
                                        Nstates = 3)
  
  # Multistate matrix: Unequal Rates
  
  Mult.Pars.Rec.UR <- asr_max_parsimony(Primates.Tree, Primates.Mult.Nost,
                                        transition_costs = Sym.Unequal.Three,
                                        Nstates = 3)
  
# Discretization of ML ancestral reconstruction values
  
  # IMPORTANT: if the difference between likelihoods at a certain node is equal #
  # or less than 0.1, both states will be considered equiprobable #
  
  DiscAceStat <- function(Node, Tolerance.Limit = 0.01){
    
    Node.Ratio <- Nodo/max(Node) 
    limit <- 1 - Tolerance.Limit/max(Node)
    Condition <- Node.Ratio >= limit
    Discrete <- paste(names(which(Condition == "TRUE")), collapse = "")
    return(Discrete)
  }
  
  # Binary matrices: Equal Rates
  
    # Rhinarium
  
    apply(Rhinarium.Anc.Rec.ER.Root, 1, DiscAceStat)
    
    # Flat nostrils
      
    apply(Flat.Nost.Anc.Rec.ER.Nodes, 1, DiscAceStat)
        
    # Downward nostrils
  
    apply(Down.Nost.Anc.Rec.ER.Root, 1, DiscAceStat)  
    
  # Binary matrices: Unequal Rates
  
    # Rhinarium
    
    apply(Rhinarium.Anc.Rec.ARD.Nodes, 1, DiscAceStat)
    
    # Flat nostrils
    
    apply(Flat.Nost.Anc.Rec.ARD.Nodes, 1, DiscAceStat) 
    
    # Downward nostrils
    
    apply(Down.Nost.Anc.Rec.ARD.Nodes, 1, DiscAceStat)
    
  # Multistate matrix: Equal Rates
  
    # Rhinarium
  
   
