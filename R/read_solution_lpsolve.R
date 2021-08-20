read_solution_lpsolve <- function(variables = variables,
                                  res = res){

  reactions <- variables$var_exp[intersect(x = which(res$solution==1),
                                           y = which(grepl(pattern =
                                                             "reaction",
                                                           x =
                                                             variables$var_exp,
                                                           fixed =
                                                             TRUE)))]

  interactions <- variables$var_exp[intersect(x =
                                                which(res$solution==1),
                                              y =
                                                which(grepl(
                                                  pattern = "interaction",
                                                  x =
                                                    variables$var_exp,
                                                  fixed =
                                                    TRUE)))]

  if(length(interactions)==0){

    stop("No solutions could be found for this setting.")

  } else {

    int <- sapply(strsplit(x = reactions,
                           split = " of ",
                           fixed = TRUE),
                  '[',
                  2)

    sif <- matrix(data = , nrow = length(interactions), ncol = 4)
    sif[, 1] <- sapply(strsplit(x =
                                  gsub(pattern =
                                         "interaction ",
                                       replacement =
                                         "",
                                       x =
                                         interactions),
                                split = "=",
                                fixed = TRUE),
                       '[',
                       1)
    sif[, 2] <- "1"
    sif[, 3] <- sapply(strsplit(x =
                                  gsub(pattern =
                                         "interaction ",
                                       replacement =
                                         "",
                                       x =
                                         interactions),
                                split = "=",
                                fixed = TRUE),
                       '[',
                       2)
    for(ii in 1:nrow(sif)){
      idx <- which(int==paste0(sif[ii, 1], "=", sif[ii, 3]))
      sif[ii, 4] <- paste0(sapply(strsplit(x = reactions[idx],
                                           split = " ",
                                           fixed = TRUE),
                                  '[',
                                  2),
                           collapse = "; ")
    }

    colnames(sif) <- c("source", "weight", "target", "reaction")

    return(sif)

  }


}
