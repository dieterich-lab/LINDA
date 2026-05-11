prepare_bn <- function(background.network = NULL,
                       as.input = NULL,
                       input.node = NULL,
                       input.scores = NULL,
                       pValThresh = NULL,
                       splice_effect_sign = NULL){

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

  targetPx<-rep(NA,2)
  sP <- get.all.shortest.paths(graph = gg,
                               from = which(rownames(adj)==input.node[1]),
                               to = which(rownames(adj)%in%tf),
                               mode = "out",
                               weights = NA)
  Px<-lapply(sP$res, function(x){return(V(gg)$name[x])})
  Px.l<-lapply(Px, length)
  Px.l<-unlist(Px.l)
  Px<-Px[which(Px.l <= 7)]
  Px<-lapply(Px, function(x){if(length(x) == 2){return(x)};
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

  if(length(Px) != 0){
    if(length(Px) > 1){
      for(i in 2:length(Px)){
        Px[[i]]<-rbind(Px[[i-1]], Px[[i]])
      }
    }
    targetPx<-rbind(targetPx, Px[[length(Px)]])
  }
  targetPx <- targetPx[2:nrow(targetPx), ]
  targetPx <- unique(targetPx)

  interactions1 <- targetPx

  ## Shortest paths between intermediate nodes
  uSpecies <- unique(c(interactions1[, 1], interactions1[, 2]))
  intNodes <- setdiff(x = uSpecies, y = c(input.node, tf))
  targetVn<-match(intNodes, V(gg)$name)

  targetPx<-rep(NA,2)
  for(ii in seq_len(length(intNodes))){

    # print(paste0("Step -- ", ii, "/", length(intNodes)))
    sP <- get.all.shortest.paths(graph = gg,
                                 from = targetVn[ii],
                                 to = targetVn[-ii],
                                 mode = "out",
                                 weights = NA)
    Px<-lapply(sP$res, function(x){return(V(gg)$name[x])})
    Px.l<-lapply(Px, length)
    Px.l<-unlist(Px.l)
    Px<-Px[which(Px.l <= 3)]
    Px<-lapply(Px, function(x){if(length(x) == 2){return(x)};
      if(length(x) == 3){return(rbind(x[1:2], x[2:3]))}
    })

    if(length(Px) != 0){
      if(length(Px) > 1){
        for(i in 2:length(Px)){
          Px[[i]]<-rbind(Px[[i-1]], Px[[i]])
        }
      }
      targetPx<-rbind(targetPx, Px[[length(Px)]])
    }
    # targetPx <- targetPx[-c(1), ]
    targetPx <- unique(targetPx)

  }

  targetPx <- unique(targetPx[2:nrow(targetPx), ])
  interactions2 <- targetPx

  interactions <- unique(rbind(interactions1, interactions2))
  reacs1 <- paste0(interactions[, 1], "=", interactions[, 2])
  reacs2 <- paste0(background.network$gene_source,
                   "=",
                   background.network$gene_target)
  idx2keep <- which(reacs2%in%reacs1)

  bg2keep <- background.network[idx2keep, ]

  ## Integrate the scores
  if(is.null(pValThresh)){

    background.network <- integrate_scores_in_bn_soft(as.input = as.input,
                                                      background.network = bg2keep,
                                                      splice_effect_sign = splice_effect_sign)

  } else {

    background.network <- integrate_scores_in_bn_hard(as.input = as.input,
                                                      background.network = bg2keep,
                                                      pValThresh = pValThresh,
                                                      splice_effect_sign = splice_effect_sign)

  }

  ## Now return the output
  returnList <- list()
  returnList[[length(returnList)+1]] <-
    background.network[complete.cases(background.network), ]
  returnList[[length(returnList)+1]] <-
    input.node
  names(returnList) <- c("background.network", "input.node")

  return(returnList)

}
