#' Performs the LINDA analysis once the inputs are given.
#'
#' @param input.scores A data-frame of two columns containing the estimated
#' transcription factor activities (TF). The data-frame should have columns ID's
#' as c('id', 'nes'), indicating to the TF ID's and the corresponding inferred
#' enrichment score.
#'
#' @param as.input a data-frame containing results from the alternative splicing/
#' exon skipping analysis. It should contain at least three columns with the
#' following ID's:
#'
#' id: a transcript/exon ID.
#'
#' effect: an exon inclusion level difference score / logFC transcript expression
#' which can be obtained from methods which are used for the detection of
#' differential AS from replicate RNA-Seq data (i.e. rMATS, MAJIQ, etc.)/methods
#' for the analysis of differential transcript expression values.
#'
#' significance: this is a significance value associated with the 'effect'
#' score. This also can be obtained from differential ES detection tools or from
#' differential transcript expression analysis.
#'
#' By default, as.input = NULL which means that LINDA will perform a
#' splice-unaware inference of upstream regulatory networks.
#'
#' @param background.network data frame of the prior knowledge network of
#' protein-protein (PPI) and domain-domain (DDI) interactions. The data-frame
#' should contain at least 6 columns with the following ID's:
#'
#' id_source: exon/transcript ID's mapping to the Pfam and protein of the source
#' interactions. In case multiple exons/transcript ID's map to a domain of a
#' specific protein, then the multiple exons/transcripts are separated by an
#' underscore symbol - '_'.
#'
#' id_target: exon/transcript ID's mapping to the Pfam and protein of the target
#' interactions. In case multiple exons/transcript ID's map to a domain of a
#' specific protein, then the multiple exons/transcripts are separated by an
#' underscore symbol - '_'.
#'
#' pfam_source: Pfam domain ID's of the interaction sources in the DDI.
#'
#' pfam_target: Pfam domain ID's of the interaction targets in the DDI.
#'
#' gene_source: gene names of the interaction sources in the PPI.
#'
#' gene_target: gene names of the interaction targets in the PPI.
#'
#' In case you are using transcripts, please load the DIGGER resource with
#' domains mapped to transcripts:
#' load(file = system.file("extdata",
#'                        "digger_human_transcripts.RData",
#'                        package = "LINDA"))
#'
#' In case you are using exons, please load the DIGGER resource with domains
#' mapped to exons:
#' load(file = system.file("extdata",
#'                        "digger_human_exons.RData",
#'                        package = "LINDA"))
#'
#' @param solverPath (optional) location path to the desired solver executable.
#' Required when solver is 'cplex' or 'cbc'.
#'
#' @param input.node (optional) a vector of perturbation inputs from where the
#' protein interactions start. By default input.node = NULL, in which case an
#' auxiliary 'Perturbation' node is added and connected on top of all uppermost
#' nodes in the background network. These are the nodes that only have outgoing
#' interactions but no incoming ones.
#'
#' @param constraints_mode character. Defines how alternative-splicing evidence
#' is integrated into the ILP. Use 'hard' to add AS constraints that force
#' significantly affected reactions to zero. Use 'soft' to keep reactions
#' feasible but penalize them in the objective function. Default is 'hard'.
#'
#' @param pValThresh (optional) a p-value/FDR threshold used in
#' constraints_mode = 'hard' to define significantly affected exons/transcripts.
#' AS-affected reactions passing this threshold and the selected
#' splice_effect_sign are excluded from the feasible ILP solution through the
#' AS constraints. By default, pValThresh = 0.05. This parameter is not used for
#' constraints_mode = 'soft'.
#'
#' @param splice_effect_sign character. Defines which sign of the AS effect is
#' considered functionally relevant. It can be 'negative', 'positive', or 'both'.
#' Default is 'negative', corresponding to exon/transcript skipping in the first
#' condition relative to the second condition when using rMATS-like
#' IncLevelDifference values.
#'
#' @param top (optional) a parameter indicating the number of TF's to be
#' considered as significantly regulated based on the absolute highest activity
#' value. The network solution must include as many of the top-regulated TF’s
#' while penalizing the inclusion of the TF’s which do not appear to be
#' regulated. By default set to top = 50. When considering all TF's, use
#' top = 'all'.
#'
#' @param lambda1 the penalization term of the primary objective of the objective
#' function - TF inclusion. This penalty factor is suggested to be set to a
#' higher value compared to other penalty parameters in order to strongly
#' penalize the inclusion of not significantly regulated TF's. By default,
#' lambda1 = 100.
#'
#' @param lambda2 the penalization term of the secondary objective of the
#' objective function - size penalty. The aim of this objective term is to
#' penalize the inclusion of spurious DDI's in the final solution. In
#' constraints_mode = 'soft', this term is additionally scaled by AS-derived
#' penalties. By default, lambda2 = 1.
#'
#' @param mipgap CPLEX parameter which sets an absolute tolerance on the gap
#' between the best integer objective and the objective of the best node
#' remaining. When this difference falls below the value of this parameter, the
#' mixed integer optimization is stopped. By default, mipgap = 0.
#'
#' @param relgap CPLEX parameter which sets a relative tolerance on the objective
#' value for the solutions in the solution pool. Solutions that are worse than
#' the incumbent solution by this measure are not kept in the solution pool. By
#' default, relgap = 0.
#'
#' @param populate CPLEX parameter which sets the maximum number of mixed integer
#' programming (MIP) solutions generated for the solution pool during each call
#' to the populate procedure. By default, populate = 500.
#'
#' @param nSolutions the number of solutions to be provided by LINDA. By default,
#' nSolutions = 100.
#'
#' @param timelimit CPLEX parameter which sets the maximum optimization time in
#' seconds. By default, timelimit = 3600.
#'
#' @param intensity CPLEX parameter which controls the trade-off between the
#' number of solutions generated for the solution pool and the amount of time or
#' memory consumed. Values from 0 to 4 invoke increasing effort to find larger
#' numbers of solutions. Higher values are more expensive in terms of time and
#' memory but are likely to yield more solutions. By default, intensity = 1.
#'
#' @param replace CPLEX parameter which designates the strategy for replacing a
#' solution in the solution pool when the solution pool has reached its capacity.
#' The value 0 replaces solutions according to a first-in, first-out policy. The
#' value 1 keeps the solutions with the best objective values. The value 2
#' replaces solutions in order to build a set of diverse solutions. By default,
#' replace = 1.
#'
#' @param threads CPLEX parameter which manages the number of threads that CPLEX
#' uses. By default, threads = 0 (let CPLEX decide).
#'
#' @param condition a parameter which can be used in the case when LINDA is
#' desired to run over multiple analyses in parallel. It is useful to distinguish
#' between the multiple ILP problems defined for each case as well as the
#' solutions obtained. By default, condition = 1.
#'
#' @param solver a character parameter indicating the type of solver to be used.
#' It can only take three different values ('lpSolve', 'cplex' and 'cbc')
#' depending on the type of solver the user wishes to use. By default,
#' solver = 'cplex'.
#'
#' @param save_res logical. If TRUE, save the result object as
#' paste0('res_', condition, '.RData'). By default, save_res = TRUE.
#'
#' @return Results list containing each of the LINDA unique solutions as well as
#' the combined ones.
#'
#' @examples
#' library(LINDA)
#' library(igraph)
#' library(XML)
#'
#' load(file = system.file("extdata", "as_data_toy.RData", package = "LINDA"))
#' load(file = system.file("extdata", "bg_toy.RData", package = "LINDA"))
#' load(file = system.file("extdata", "tf_act_toy.RData", package = "LINDA"))
#'
#' \dontrun{
#' # Default: hard AS constraints.
#' res_hard <- runLINDA(input.scores = input.scores,
#'                      as.input = as.input,
#'                      background.network = bg,
#'                      solverPath = "~/Downloads/cplex",
#'                      input.node = NULL,
#'                      constraints_mode = "hard",
#'                      pValThresh = 0.05,
#'                      splice_effect_sign = "negative",
#'                      top = 2,
#'                      lambda1 = 10,
#'                      lambda2 = 0.001,
#'                      mipgap = 0.001,
#'                      relgap = 0.001)
#'
#' # Soft AS penalties instead of hard AS constraints.
#' res_soft <- runLINDA(input.scores = input.scores,
#'                      as.input = as.input,
#'                      background.network = bg,
#'                      solverPath = "~/Downloads/cplex",
#'                      constraints_mode = "soft",
#'                      splice_effect_sign = "negative",
#'                      top = 2,
#'                      lambda1 = 10,
#'                      lambda2 = 0.001)
#' }
#'
#' @export

runLINDA <- function(input.scores = input.scores,
                     as.input = NULL,
                     background.network = bg,
                     solverPath = NULL,
                     input.node = NULL,
                     constraints_mode = "hard",
                     pValThresh = 0.05,
                     splice_effect_sign = "negative",
                     top = 50,
                     lambda1 = 100,
                     lambda2 = 1,
                     mipgap = 0,
                     relgap = 0,
                     populate = 500,
                     nSolutions = 100,
                     timelimit = 3600,
                     intensity = 1,
                     replace = 1,
                     threads = 0,
                     condition = 1,
                     solver = "cplex",
                     save_res = TRUE){

  options(scipen = 999)

  constraints_mode <- match.arg(constraints_mode, choices = c("hard", "soft"))
  splice_effect_sign <- match.arg(splice_effect_sign,
                                  choices = c("negative", "positive", "both"))

  checkInputs(input.scores = input.scores,
              as.input = as.input,
              background.network = background.network,
              solverPath = solverPath,
              input.node = input.node,
              constraints_mode = constraints_mode,
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
              solver = solver,
              save_res = save_res)

  input.scores <- bin_measurements(input.scores = input.scores, top = top)

  colnames(background.network)[which(colnames(background.network) == "id_source")] <- "exon_source"
  colnames(background.network)[which(colnames(background.network) == "id_target")] <- "exon_target"

  bn <- prepare_bn(background.network = background.network,
                   as.input = as.input,
                   input.node = input.node,
                   input.scores = input.scores,
                   constraints_mode = constraints_mode,
                   pValThresh = pValThresh,
                   splice_effect_sign = splice_effect_sign)

  variables <- create_variables(background.network = bn$background.network)

  if(constraints_mode == "soft"){

    objective.function <- write_objective_function_soft(background.network =
                                                          bn$background.network,
                                                        variables = variables,
                                                        input.scores = input.scores,
                                                        lambda1 = lambda1,
                                                        lambda2 = lambda2)

  } else {

    objective.function <- write_objective_function(background.network =
                                                     bn$background.network,
                                                   variables = variables,
                                                   input.scores = input.scores,
                                                   lambda1 = lambda1,
                                                   lambda2 = lambda2)

  }

  allC <- write_all_constraints(variables = variables,
                                bn = bn,
                                constraints_mode = constraints_mode,
                                pValThresh = pValThresh,
                                input.scores = input.scores)

  if(solver == "cplex"){

    res <- computeILP_CPLEX(variables = variables,
                            background.network = bn$background.network,
                            solverPath = solverPath,
                            mipgap = mipgap,
                            relgap = relgap,
                            populate = populate,
                            nSolutions = nSolutions,
                            timelimit = timelimit,
                            intensity = intensity,
                            replace = replace,
                            condition = condition,
                            threads = threads,
                            allC = allC,
                            objective.function = objective.function)

  } else {

    if(solver == "lpSolve"){

      res <- computeILP_LpSolve(variables = variables,
                                objective.function = objective.function,
                                allC = allC,
                                condition = condition)

    } else {

      res <- computeILP_CBC(variables = variables,
                            objective.function = objective.function,
                            background.network = bn$background.network,
                            solverPath = solverPath,
                            mipgap = mipgap,
                            relgap = relgap,
                            populate = populate,
                            nSolutions = nSolutions,
                            timelimit = timelimit,
                            intensity = intensity,
                            replace = replace,
                            condition = condition,
                            threads = threads,
                            allC = allC)

    }

  }

  if(save_res){
    save(res, file = paste0("res_", condition, ".RData"))
  }

  return(res)

}
