read_solution_cplex <- function(variables = variables,
                                background.network = background.network,
                                condition = 1){

  cplexSolutionFileName <- paste0("results_", condition, ".txt")

  reacVar <- variables$var[which(grepl(pattern = "reaction",
                                       x = variables$var_exp))]
  intVar <- variables$var[which(grepl(pattern = "interaction",
                                      x = variables$var_exp))]

  cplexSolutionData <- xmlParse(cplexSolutionFileName)
  cplexSolution <- xmlToList(cplexSolutionData)

  cplexSolutionIntAll <- list()
  cplexSolutionReacAll <- list()

  sifAll <- list()

  # 2:(length(cplexSolution)-1)
  for(ii in 2:(length(cplexSolution)-1)){

    currSolution <- cplexSolution[[ii]][[4]]

    varvar <- unlist(lapply(currSolution, '[', 1))
    valval <- as.numeric(unlist(lapply(currSolution, '[', 3)))

    idxReac <- intersect(x = which(varvar%in%reacVar), y = which(valval>=0.99))
    # idxInt <- intersect(x = which(varvar%in%intVar), y = which(valval>=0.99))

    if(length(idxReac)>0){

      reactions <-
        sapply(strsplit(x =
                          variables$var_exp[
                            which(variables$var%in%varvar[idxReac])],
                        split = " ", fixed = TRUE), '[', 2)
      interactions <-
        sapply(strsplit(x =
                          variables$var_exp[
                            which(variables$var%in%varvar[idxReac])],
                        split = " ", fixed = TRUE), '[', 4)

      uInt <- unique(interactions)
      currSIF <- matrix(data = , nrow = length(uInt), ncol = 4)
      currSIF[, 1] <- sapply(strsplit(x = uInt,
                                      split = "=",
                                      fixed = TRUE),
                             '[',
                             1)
      currSIF[, 2] <- "1"
      currSIF[, 3] <- sapply(strsplit(x = uInt,
                                      split = "=",
                                      fixed = TRUE),
                             '[',
                             2)
      for(jj in seq_len(length(uInt))){

        idx <-
          which(grepl(pattern =
                        paste0(" ",
                               uInt[jj]),
                      x = variables$var_exp[
                        which(variables$var%in%varvar[idxReac])],
                      fixed = TRUE))
        currSIF[jj, 4] <- paste0(reactions[idx], collapse = "; ")

      }

      sifAll[[length(sifAll)+1]] <- currSIF

    }

  }

  for(ii in seq_len(length(sifAll))){

    if(ii==1){

      combSIF <- sifAll[[ii]][, 1:3]

    } else {

      currSIF <- sifAll[[ii]][, 1:3]
      for(jj in seq_len(nrow(currSIF))){

        idx1 <- which(combSIF[, 1]==currSIF[jj, 1])
        idx2 <- which(combSIF[, 3]==currSIF[jj, 3])
        idx <- intersect(x = idx1, y = idx2)
        if(length(idx)>0){

          combSIF[idx, 2] <- as.character(as.numeric(combSIF[idx, 2])+1)

        } else {

          combSIF <- rbind(combSIF, t(as.matrix(currSIF[jj, ])))

        }

      }

    }

  }

  domains <- rep("", nrow(combSIF))
  for(ii in seq_len(nrow(combSIF))){

    ss <- combSIF[ii, 1]
    tt <- combSIF[ii, 3]

    ud <- c()
    for(jj in seq_len(length(sifAll))){

      idx1 <- which(sifAll[[jj]][, 1]==ss)
      idx2 <- which(sifAll[[jj]][, 3]==tt)
      idx <- intersect(x = idx1, y = idx2)
      if(length(idx)>0){

        ud <- c(ud, strsplit(x = sifAll[[jj]][idx, 4],
                             split = "; ",
                             fixed = TRUE)[[1]])

      }

    }

    domains[ii] <- paste0(unique(ud), collapse = "; ")
  }

  sif <- matrix(data = , nrow = nrow(combSIF), ncol = 4)
  sif[, 1:3] <- combSIF
  sif[, 4] <- domains

  idx2rem <- setdiff(x = which(sif[, 4]==""),
                     y = which(sif[, 1]=="Perturbation"))
  if(length(idx2rem)>0){
    sif <- sif[-idx2rem, ]
  }

  if(nrow(sif)>0){

    colnames(sif) <- c("source", "weight", "target", "reaction")

    # return(sif)
    returnList <- list()
    returnList[[length(returnList)+1]] <- sifAll
    returnList[[length(returnList)+1]] <- sif

    names(returnList) <- c("separate_interactions", "combined_interactions")

    return(returnList)
  } else {

    return(NULL)

  }

}
