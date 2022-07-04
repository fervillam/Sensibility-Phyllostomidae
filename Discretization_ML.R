function(Reconstruction.Matrix, Tolerance.Limit = 0.01){
  Nodes <- dim(Reconstruction.Matrix)[1]
  States <- dim(Reconstruction.Matrix)[2]
  
  Limit <- Tolerance.Limit * 100
  
  Lik.Limit <- list()
  
  for(i in 1:Nodes){
    out <- (Reconstruction.Matrix[i, ])/max(Reconstruction.Matrix[i, ])
    Lik.Limit[[i]] <- out
  }

  if(Lik.Limit[1:Nodes] >= Limit, "TRUE"){
    State.Vector <- ()
  }
  
