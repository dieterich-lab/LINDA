computeILP_CBC <- function(variables = variables,
                           objective.function = objective.function,
                           background.network = background.network,
                           solverPath = solverPath,
                           mipgap = 0.05,
                           relgap = 0.05,
                           populate = 1000,
                           nSolutions = 100,
                           timelimit = 3600,
                           intensity = 2,
                           replace = 1,
                           condition = condition,
                           threads = threads,
                           allC = allC){

  bounds <- write_bounds(variables = variables)
  binaries <- write_binaries(variables = variables)

  # write the .lp file
  data = paste0("testFile_", condition, ".lp")
  write("enter Problem", data)
  write("", data, append = TRUE)
  write("Minimize", data, append = TRUE)
  write(objective.function, data, append = TRUE)
  write("Subject To", data, append = TRUE)
  write(allC, data, append = TRUE)
  write("Bounds", data, append = TRUE)
  write(bounds, data, append = TRUE)
  write("Binaries", data, append = TRUE)
  write(binaries, data, append = TRUE)
  write("End", data, append = TRUE)

  solve_with_cbc(solverPath = solverPath, condition = condition,
                 timelimit = timelimit, relgap = relgap)
  
  # Read, cleanup and return solution
  sif <- read_solution_cbc(variables = variables, condition = condition)
  
  cleanupILP(condition = condition)
  
  return(sif)

}
