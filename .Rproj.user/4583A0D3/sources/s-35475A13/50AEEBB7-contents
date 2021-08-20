solve_with_cbc <- function(solverPath = solverPath,
                           condition = 1,
                           timelimit = timelimit,
                           relgap = relgap){

  cbc_command <- paste0(solverPath, " testFile_", condition,
                        ".lp -seconds ", timelimit,
                        " -ratio ", relgap,
                        " solve printi csv solu ",
                        paste0("solution_", condition, ".csv"))

  system(cbc_command)

}
