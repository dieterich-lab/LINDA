checkInputs <- function(input.scores = input.scores,
                        as.input = as.input,
                        background.network = bg,
                        solverPath = solverPath,
                        input.node = input.node,
                        constraints_mode = "hard",
                        pValThresh = 0.05,
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
                        solver = "lpSolve",
                        save_res = TRUE){

  # input.scores
  if(!is.data.frame(input.scores) ||
     ncol(input.scores) != 2 ||
     length(intersect(x = colnames(input.scores),
                      y = c("id", "nes"))) < 2){

    stop("The 'input.scores' object should be a data-frame with two columns and
         have c('id', 'nes') as column ID's. Please check your inputs.")

  }

  # as.input
  if(!is.null(as.input)){

    if(!is.data.frame(as.input) ||
       ncol(as.input) < 3 ||
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
     ncol(background.network) < 6 ||
     length(intersect(x = colnames(background.network),
                      y = c("gene_source", "gene_target",
                            "pfam_source", "pfam_target",
                            "id_source", "id_target"))) < 6){

    stop("The 'background.network' object should be a data-frame with at least 6
         columns and have at least c('gene_source', 'gene_target',
         'pfam_source', 'pfam_target', 'id_source', 'id_target') as column
         ID's. Please check your inputs.")

  }

  # input.node
  if(!is.null(input.node)){

    if(!is.character(input.node)){

      stop("the 'input.node' object should either be NULL or a vector of
         characters. Please check your inputs.")

    }

  }

  # constraints mode
  if(length(constraints_mode) != 1 ||
     !is.character(constraints_mode) ||
     !(constraints_mode %in% c("hard", "soft"))){

    stop("The 'constraints_mode' parameter should be either 'hard' or 'soft'.
         Please check your inputs.")

  }

  # pValThresh parameter
  if(!is.null(pValThresh)){

    if(!is.numeric(pValThresh) ||
       length(pValThresh) != 1 ||
       is.na(pValThresh) ||
       pValThresh < 0 ||
       pValThresh > 1){

      stop("The 'pValThresh' parameter should be NULL or a numeric value between
           0 and 1. Please check your inputs.")

    }

  }

  if(constraints_mode == "hard" && !is.null(as.input) && is.null(pValThresh)){

    stop("When constraints_mode = 'hard' and as.input is provided, 'pValThresh'
         should be a numeric value between 0 and 1. Please check your inputs.")

  }

  # solver
  if(length(solver) != 1 ||
     !is.character(solver) ||
     !(solver %in% c("lpSolve", "cplex", "cbc"))){

    stop("The 'solver' parameter should be a character of 'lpSolve', 'cplex' or
         'cbc' value. Please check your inputs.")

  }

  # splice_effect_sign parameter
  if(length(splice_effect_sign) != 1 ||
     !is.character(splice_effect_sign) ||
     !(splice_effect_sign %in% c("positive", "negative", "both"))){

    stop("The 'splice_effect_sign' parameter should be a character of 'positive',
         'negative' or 'both' value. Please check your inputs.")

  }

  if(!(identical(top, "all") ||
       (is.numeric(top) && length(top) == 1 && !is.na(top) &&
        top == round(top) && top > 0))){

    stop("The 'top' parameter should either be a positive numeric integer or
         'all'. Please check your inputs.")

  }

  if(is.numeric(top) && top > nrow(input.scores)){

    stop("The 'top' parameter cannot be larger than the number of rows in
         'input.scores'. Please check your inputs.")

  }

  if(!is.numeric(lambda1) || length(lambda1) != 1 || is.na(lambda1)){

    stop("The 'lambda1' parameter should be numeric. Please check your inputs.")

  }

  if(!is.numeric(lambda2) || length(lambda2) != 1 || is.na(lambda2)){

    stop("The 'lambda2' parameter should be numeric. Please check your inputs.")

  }

  if(!is.numeric(mipgap) || length(mipgap) != 1 || is.na(mipgap) ||
     mipgap > 1 ||
     mipgap < 0){

    stop("The 'mipgap' parameter should be numeric between 0 and 1. Please check
         your inputs.")

  }

  if(!is.numeric(relgap) || length(relgap) != 1 || is.na(relgap) ||
     relgap > 1 ||
     relgap < 0){

    stop("The 'relgap' parameter should be numeric between 0 and 1. Please check
         your inputs.")

  }

  if(!is.numeric(populate) || length(populate) != 1 || is.na(populate) ||
     populate != round(populate) ||
     populate <= 0){

    stop("The 'populate' parameter should be a positive numeric integer. Please
         check your inputs.")

  }

  if(!is.numeric(nSolutions) || length(nSolutions) != 1 || is.na(nSolutions) ||
     nSolutions != round(nSolutions) ||
     nSolutions <= 0){

    stop("The 'nSolutions' parameter should be a positive numeric integer.
         Please check your inputs.")

  }

  if(!is.numeric(timelimit) || length(timelimit) != 1 || is.na(timelimit) ||
     timelimit <= 0){

    stop("The 'timelimit' parameter should be a positive numeric value. Please
         check your inputs.")

  }

  if(length(intensity) != 1 || is.na(intensity) ||
     length(intersect(x = intensity, y = 0:4)) != 1){

    stop("The 'intensity' parameter should be numeric value between 0 and 4.
         Please check your inputs.")

  }

  if(length(replace) != 1 || is.na(replace) ||
     length(intersect(x = replace, y = 0:2)) != 1){

    stop("The 'replace' parameter should be numeric value between 0 and 2.
         Please check your inputs.")

  }

  if(!is.numeric(threads) || length(threads) != 1 || is.na(threads) ||
     threads != round(threads) ||
     threads < 0){

    stop("The 'threads' parameter should be a positive (or 0) numeric integer.
         Please check your inputs.")

  }

  if(!is.numeric(condition) || length(condition) != 1 || is.na(condition) ||
     condition != round(condition) ||
     condition <= 0){

    stop("The 'condition' parameter should be a positive numeric integer. Please
         check your inputs.")

  }

  if(!is.logical(save_res) || length(save_res) != 1 || is.na(save_res)){

    stop("The 'save_res' parameter should be TRUE or FALSE. Please check your
         inputs.")

  }

  # solverPath
  if(solver != "lpSolve"){

    if(is.null(solverPath)){

      stop(paste0("You are trying to use ", solver, " as a solver. Please
                  provide a path pointing to the ", solver, " executable solver
                  file in the 'solverPath' parameter."))

    } else {

      if(!file.exists(solverPath)){

        stop("The path to the solver that you provided seems to not exist.
             Please check your inputs.")

      }

    }

  }

}
