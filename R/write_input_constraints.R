write_input_constraints <- function(variables = variables,
                                    input.node = input.node){

  # print("Writing Constraints 1/6..")
  idx <- which(variables$var_exp%in%paste0("node ", input.node))

  cc <- rep("", length(idx))
  for(ii in 1:length(idx)){

    cc[ii] <- paste0(variables$var[idx[ii]], " = 1")

  }

  return(cc)

}
