prepare_bn <- function(background.network = NULL,
                       as.input = NULL,
                       input.node = NULL,
                       input.scores = NULL){

  print("Processing the Background Network..")

  # remove first self-interactions
  idx2rem <-
    which(background.network$gene_source==background.network$gene_target)
  if(length(idx2rem)>0){
    background.network <- background.network[-idx2rem, ]
  }

  if(is.null(input.node)){

    ppi <- unique(background.network[, c("gene_source", "gene_target")])

    top.proteins <- unique(setdiff(x = ppi$gene_source, y = ppi$gene_target))

    for(ii in seq_len(length(top.proteins))){

      upfam <-
        unique(background.network$pfam_source[
          which(background.network$gene_source==top.proteins[ii])])
      utrans <-
        unique(background.network$exon_source[
          which(background.network$gene_source==top.proteins[ii])])

      toBind <- matrix(data = ,
                       nrow = length(upfam),
                       ncol = ncol(background.network))
      colnames(toBind) <- colnames(background.network)
      for(jj in seq_len(nrow(toBind))){

        toBind[jj, 1] <- "Perturbation"
        toBind[jj, 2] <- utrans[jj]
        toBind[jj, 3] <- "Perturbation"
        toBind[jj, 4] <- upfam[jj]
        toBind[jj, 5] <- "Perturbation"
        toBind[jj, 6] <- top.proteins[ii]

      }

      toBind <- as.data.frame(toBind)
      background.network <- rbind(background.network, toBind)

    }

    background.network <- unique(background.network)
    input.node <- "Perturbation"

  }

  ## Shortest paths from the inputs to the measurements
  df <- unique(background.network[, c("gene_source", "gene_target")])
  gg <- graph_from_data_frame(d = df, directed = TRUE)
  adj <- get.adjacency(graph = gg)
  tf <- input.scores$id

  targetPaths<-rep(NA,2)
  sP <- get.all.shortest.paths(graph = gg,
                               from = which(rownames(adj)==input.node[1]),
                               to = which(rownames(adj)%in%tf),
                               mode = "out",
                               weights = NA)
  Paths<-lapply(sP$res, function(x){return(V(gg)$name[x])})
  Paths.l<-lapply(Paths, length)
  Paths.l<-unlist(Paths.l)
  Paths<-Paths[which(Paths.l <= 7)]
  Paths<-lapply(Paths, function(x){if(length(x) == 2){return(x)};
    if(length(x) == 3){return(rbind(x[1:2], x[2:3]))};
    if(length(x) == 4){
      return(rbind(x[1:2], x[2:3], x[3:4]))
    };
    if(length(x) == 6){
      return(rbind(x[1:2], x[2:3], x[3:4], x[4:5], x[5:6]))
    };
    if(length(x) == 7){
      return(rbind(x[1:2], x[2:3], x[3:4], x[4:5], x[5:6], x[6:7]))
    };
    if(length(x) == 5){
      return(rbind(x[1:2], x[2:3], x[3:4], x[4:5]))}
  })

  if(length(Paths) != 0){
    if(length(Paths) > 1){
      for(i in 2:length(Paths)){
        Paths[[i]]<-rbind(Paths[[i-1]], Paths[[i]])
      }
    }
    targetPaths<-rbind(targetPaths, Paths[[length(Paths)]])
  }
  targetPaths <- targetPaths[2:nrow(targetPaths), ]
  targetPaths <- unique(targetPaths)

  interactions1 <- targetPaths

  ## Shortest paths between intermediate nodes
  uSpecies <- unique(c(interactions1[, 1], interactions1[, 2]))
  intNodes <- setdiff(x = uSpecies, y = c(input.node, tf))
  targetVn<-match(intNodes, V(gg)$name)

  targetPaths<-rep(NA,2)
  for(ii in seq_len(length(intNodes))){

    # print(paste0("Step -- ", ii, "/", length(intNodes)))
    sP <- get.all.shortest.paths(graph = gg,
                                 from = targetVn[ii],
                                 to = targetVn[-ii],
                                 mode = "out",
                                 weights = NA)
    Paths<-lapply(sP$res, function(x){return(V(gg)$name[x])})
    Paths.l<-lapply(Paths, length)
    Paths.l<-unlist(Paths.l)
    Paths<-Paths[which(Paths.l <= 3)]
    Paths<-lapply(Paths, function(x){if(length(x) == 2){return(x)};
      if(length(x) == 3){return(rbind(x[1:2], x[2:3]))}
    })

    if(length(Paths) != 0){
      if(length(Paths) > 1){
        for(i in 2:length(Paths)){
          Paths[[i]]<-rbind(Paths[[i-1]], Paths[[i]])
        }
      }
      targetPaths<-rbind(targetPaths, Paths[[length(Paths)]])
    }
    # targetPaths <- targetPaths[-c(1), ]
    targetPaths <- unique(targetPaths)

  }

  targetPaths <- unique(targetPaths[2:nrow(targetPaths), ])
  interactions2 <- targetPaths

  interactions <- unique(rbind(interactions1, interactions2))
  reacs1 <- paste0(interactions[, 1], "=", interactions[, 2])
  reacs2 <- paste0(background.network$gene_source,
                   "=",
                   background.network$gene_target)
  idx2keep <- which(reacs2%in%reacs1)

  bg2keep <- background.network[idx2keep, ]

  ## Integrate the scores
  background.network <- integrate_scores_in_bn(as.input = as.input,
                                               background.network = bg2keep)

  ## Now return the output
  returnList <- list()
  returnList[[length(returnList)+1]] <-
    background.network[complete.cases(background.network), ]
  returnList[[length(returnList)+1]] <-
    input.node
  names(returnList) <- c("background.network", "input.node")

  return(returnList)

}
