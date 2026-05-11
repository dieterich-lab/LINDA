checkInputs <- function(input.scores = input.scores,
                        as.input = as.input,
                        background.network = bg,
                        solverPath = solverPath,
                        input.node = input.node,
                        pValThresh = pValThresh,
                        splice_effect_sign = splice_effect_sign,
                        top = top,
                        lambda1 = lambda1,
                        lambda2 = lambda2,
                        mipgap = mipgap,
                        relgap = relgap,
                        populate = populate,
                        nSolutions = nSolutions,
                        timelimit = timelimit,
                        intensity = intensity,
                        replace = replace,
                        threads = threads,
                        condition = condition,
                        solver = "lpSolve"){

  # input.scores
  if(!is.data.frame(input.scores) ||
     ncol(input.scores)!=2 ||
     length(intersect(x = colnames(input.scores),
                      y = c("id", "nes"))) < 2){

    stop("The 'input.scores' object should be a data-frame with two columns and
         have c('id', 'nes') as column ID's. Please check your inputs.")

  }

  # as.input
  if(!is.null(as.input)){

    if(!is.data.frame(as.input) ||
       ncol(as.input)<3 ||
       length(intersect(x = colnames(as.input),
                        y = c("id",
                              "effect",
                              "significance"))) < 3){

      stop("The 'as.input' object should either be NULL or a data-frame with
            three columns and have at least c('id','effect', 'significance')
            as column ID's. Please check your inputs.")

    }

  }

  # background.network
  if(!is.data.frame(background.network) ||
     ncol(background.network)<6 ||
     length(intersect(x = colnames(background.network),
                      y = c("gene_source", "gene_target",
                            "pfam_source", "pfam_target",
                            "id_source", "id_target"))) < 6){

    stop("The 'background.network' object should be a data-frame with at least 6
         columns and have at lease c('gene_source', 'gene_target',
         'pfam_source', 'pfam_target', 'exon_source', 'exon_target') as column
         ID's. Please check your inputs.")

  }

  # input.node
  if(!is.null(input.node)){

    if(!is.character(input.node)){

      stop("the 'input.node' object should either be NULL or a vector of
         characters. Please check your inputs.")

    }

  }

  # CPLEX parameters
  if(length(intersect(x = solver, y = c("lpSolve", "cplex", "cbc"))) != 1 ||
     !is.character(solver)){

    stop("The 'solver' parameter should be a character of 'lpSolve', 'cplex' or
         'cbc' value. Please check your inputs.")

  }

  # if(!is.numeric(pValThresh) ||
  #    length(pValThresh) != 1 || !is.null(pValThresh)){
  #
  #   stop("The 'pValThresh' parameter should either be numeric for the
  #   hard-constrained version of LINDA or set to NULL for the soft-constrained
  #   analyses. Please check your inputs.")
  #
  # }

  # splice_effect_sign parameters
  if(length(intersect(x = splice_effect_sign, y = c("positive", "negative", "both"))) != 1 ||
     !is.character(splice_effect_sign)){

    stop("The 'splice_effect_sign' parameter should be a character of 'positive',
         'negative' or 'both' value. Please check your inputs.")

  }

  if(!is.numeric(top)){

    stop("The 'top' parameter should be numeric. Please check your inputs.")

  }

  if(!is.numeric(lambda1)){

    stop("The 'lambda1' parameter should be numeric. Please check your inputs.")

  }

  if(!is.numeric(lambda2)){

    stop("The 'lambda2' parameter should be numeric. Please check your inputs.")

  }

  if(!is.numeric(mipgap) ||
     mipgap > 1 ||
     mipgap < 0){

    stop("The 'mipgap' parameter should be numeric between 0 and 1. Please check
         your inputs.")

  }

  if(!is.numeric(relgap) ||
     relgap > 1 ||
     relgap < 0){

    stop("The 'relgap' parameter should be numeric between 0 and 1. Please check
         your inputs.")

  }

  if(!is.numeric(populate) ||
     populate!=round(populate) ||
     populate <= 0){

    stop("The 'populate' parameter should be a positive numeric integer. Please
         check your inputs.")

  }

  if(!is.numeric(nSolutions) ||
     nSolutions!=round(nSolutions) ||
     nSolutions <= 0){

    stop("The 'nSolutions' parameter should be a positive numeric integer.
         Please check your inputs.")

  }

  if(!is.numeric(timelimit) ||
     timelimit<=0){

    stop("The 'timelimit' parameter should be a positive numeric value. Please
         check your inputs.")

  }

  if(length(intersect(x = intensity, y = 0:4)) != 1){

    stop("The 'intensity' parameter should be numeric value between 0 and 4.
         Please check your inputs.")

  }

  if(length(intersect(x = replace, y = 0:2)) != 1){

    stop("The 'replace' parameter should be numeric value between 0 and 2.
         Please check your inputs.")

  }

  if(!is.numeric(threads) ||
     threads!=round(threads) ||
     threads < 0){

    stop("The 'nSolutions' parameter should be a positive (or 0) numeric
         integer. Please check your inputs.")

  }

  if(!is.numeric(condition) ||
     condition!=round(condition) ||
     condition <= 0){

    stop("The 'condition' parameter should be a positive numeric integer. Please
         check your inputs.")

  }

  # solverPath
  if(solver != "lpSolve"){

    if(is.null(solverPath)){

      stop(paste0("You are trying to use ", solver, " as a solver. So please
                  insert a path pointing to the ", solver, " executable solver
                  file in the 'solver' parameter."))

    } else {

      if(!file.exists(solverPath)){

        stop("The path to the solver that you provided seem to not exist. Please
             check your inputs.")

      }

    }

  }

}
