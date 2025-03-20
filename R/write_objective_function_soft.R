write_objective_function_soft <- function(background.network = background.network,
                                          variables = variables,
                                          input.scores = input.scores,
                                          lambda1=10,
                                          lambda2=0.01){
  
  print("Writing the objective function and constraints.
        This might take a bit of time..")
  
  # Write main objective - Minimizing the TF scores
  tfNodes <- input.scores$id
  idx <- which(variables$var_exp%in%paste0("node ", tfNodes))
  if(length(idx)==0){
    stop("No input protein mapped to the background network.
         Check inputs/background network")
  }
  
  objective.function <- ""
  cnt <- 1
  for(ii in seq_len(nrow(input.scores))){
    
    idx <- which(variables$var_exp==paste0("node ", input.scores$id[ii]))
    if(length(idx)==1){
      
      if(cnt == 1){
        
        if(input.scores$bin[ii]==1){
          objective.function <- paste0(objective.function,
                                       lambda1,
                                       " ",
                                       variables$var[idx])
        } else {
          objective.function <- paste0(objective.function,
                                       "- ",
                                       lambda1,
                                       " ",
                                       variables$var[idx])
        }
        cnt <- cnt + 1
        
      } else {
        
        if(input.scores$bin[ii]==1){
          objective.function <- paste0(objective.function,
                                       " + ",
                                       lambda1,
                                       " ",
                                       variables$var[idx])
        } else {
          objective.function <- paste0(objective.function,
                                       " - ",
                                       lambda1,
                                       " ",
                                       variables$var[idx])
        }
        
      }
      
    }
    
  }
  
  # Write third objective - size penalty factor
  idx <- which(grepl(pattern = "reaction", x = variables$var_exp, fixed = TRUE))
  if(lambda2>0){
    obj <- paste0(" + ", lambda2*background.network$min_score + lambda2, " ", variables$var[idx])
    obj <- paste0(obj, collapse = "")
    objective.function <- paste0(objective.function, obj)
  } else {
    if(lambda2<0){
      obj <- paste0(" - ", (abs(lambda2)*background.network$min_score + abs(lambda2)), " ", variables$var[idx])
      obj <- paste0(obj, collapse = "")
      objective.function <- paste0(objective.function, obj)
    }
  }
  
  return(objective.function)
  
}
