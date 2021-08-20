write_all_constraints <- function(variables = variables,
                                  bn = bn,
                                  pValThresh = pValThresh,
                                  input.scores = input.scores){

  cc1 <- write_input_constraints(variables = variables,
                                 input.node = bn$input.node)

  cc2 <- write_as_constraints(background.network = bn$background.network,
                              variables = variables,
                              pValThresh = pValThresh)

  cc3 <- write_mediation_constraints(variables = variables,
                                     background.network = bn$background.network)

  cc4 <- write_inout_constraints(background.network = bn$background.network,
                                 input.node = bn$input.node,
                                 input.scores = input.scores,
                                 variables = variables)

  cc5  <- write_domain_constraints(background.network = bn$background.network,
                                   variables = variables)

  cc6 <- write_loop_constraints(variables = variables,
                                background.network = bn$background.network)

  allC <- c(cc1, cc2, cc3, cc4, cc5, cc6)
  idx2rem <- c(which(allC==""), which(is.na(allC)))
  if(length(idx2rem)>0){
    allC <- allC[-idx2rem]
  }

  return(allC)

}
