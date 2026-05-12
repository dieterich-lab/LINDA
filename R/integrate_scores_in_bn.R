integrate_scores_in_bn <- function(as.input = as.input,
                                   background.network = background.network,
                                   pValThresh = NULL,
                                   splice_effect_sign = splice_effect_sign,
                                   constraints_mode = ifelse(is.null(pValThresh), "soft", "hard")){

  constraints_mode <- match.arg(constraints_mode, choices = c("hard", "soft"))

  if(constraints_mode == "soft"){

    return(integrate_scores_in_bn_soft(as.input = as.input,
                                       background.network = background.network,
                                       splice_effect_sign = splice_effect_sign))

  }

  return(integrate_scores_in_bn_hard(as.input = as.input,
                                     background.network = background.network,
                                     pValThresh = pValThresh,
                                     splice_effect_sign = splice_effect_sign))

}
