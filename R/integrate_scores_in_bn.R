integrate_scores_in_bn <- function(as.input = as.input,
                                   background.network = background.network,
                                   pValThresh = pValThresh,
                                   splice_effect_sign = splice_effect_sign){

  if(!is.null(as.input)){
    colnames(as.input) <- c("exon_id", "IncLevelDifference", "pval")
  }

  # print("Integrating AS scores in the Background Network..")
  source_score <- rep(0, nrow(background.network))
  target_score <- rep(0, nrow(background.network))
  source_fdr <- rep(1, nrow(background.network))
  target_fdr <- rep(1, nrow(background.network))

  if(is.null(as.input)){

    min_score <- pmin(source_score, target_score)
    min_fdr <- pmin(source_fdr, target_fdr)

    background.network$source_score <- source_score
    background.network$target_score <- target_score
    background.network$min_score <- min_score

    background.network$source_fdr <- source_fdr
    background.network$target_fdr <- target_fdr
    background.network$min_fdr <- min_fdr

  } else {

    as <- as.input

    uTranscripts <- unique(as$exon_id)
    for(ii in seq_len(length(uTranscripts))){

      # fdr <- min(as$pval[which(as$exon_id==uTranscripts[ii])],
      #            na.rm = TRUE)
      # score <- min(as$IncLevelDifference[intersect(
      #   x = which(as$exon_id==uTranscripts[ii]),
      #   y = which(as$pval==fdr))],
      #   na.rm = TRUE)

      fdr <- fisher(as$pval[which(as$exon_id==uTranscripts[ii])])
      score <- mean(as$IncLevelDifference[intersect(
        x = which(as$exon_id==uTranscripts[ii]),
        y = which(as$pval==fdr))],
        na.rm = TRUE)

      idx <- which(grepl(pattern = uTranscripts[ii],
                         x = background.network$exon_source,
                         fixed = TRUE))
      if(length(idx)>0){
        for(jj in seq_len(length(idx))){
          if(fdr < source_fdr[idx[jj]]){
            source_score[idx[jj]] <- score
            source_fdr[idx[jj]] <- fdr
          }
        }
      }

      idx <- which(grepl(pattern = uTranscripts[ii],
                         x = background.network$exon_target,
                         fixed = TRUE))
      if(length(idx)>0){
        for(jj in seq_len(length(idx))){
          if(fdr < target_fdr[idx[jj]]){
            target_score[idx[jj]] <- score
            target_fdr[idx[jj]] <- fdr
          }
        }
      }

    }

    if(splice_effect_sign == "negative"){

      min_score <- rep("", nrow(background.network))
      min_fdr <- rep("", nrow(background.network))
      for(ii in 1:nrow(background.network)){

        if(source_score[ii]<0 && source_fdr[ii]<pValThresh){
          min_score[ii] <- -1
          min_fdr[ii] <- source_fdr[ii]
        } else {
          if(target_score[ii]<0 && target_fdr[ii]<pValThresh){
            min_score[ii] <- -1
            min_fdr[ii] <- target_fdr[ii]
          } else {
            min_score[ii] <- 1
            min_fdr[ii] <- max(source_fdr[ii], target_fdr[ii])
          }
        }
      }

    } else {

      if(splice_effect_sign == "positive"){

        min_score <- rep("", nrow(background.network))
        min_fdr <- rep("", nrow(background.network))
        for(ii in 1:nrow(background.network)){

          if(source_score[ii]>0 && source_fdr[ii]<pValThresh){
            min_score[ii] <- -1
            min_fdr[ii] <- source_fdr[ii]
          } else {
            if(target_score[ii]>0 && target_fdr[ii]<pValThresh){
              min_score[ii] <- -1
              min_fdr[ii] <- target_fdr[ii]
            } else {
              min_score[ii] <- 1
              min_fdr[ii] <- max(source_fdr[ii], target_fdr[ii])
            }
          }
        }

      } else {

        if(splice_effect_sign == "both"){

          min_score <- rep("", nrow(background.network))
          min_fdr <- rep("", nrow(background.network))
          for(ii in 1:nrow(background.network)){

            if(source_fdr[ii]<pValThresh){
              min_score[ii] <- -1
              min_fdr[ii] <- source_fdr[ii]
            } else {
              if(target_fdr[ii]<pValThresh){
                min_score[ii] <- -1
                min_fdr[ii] <- target_fdr[ii]
              } else {
                min_score[ii] <- 1
                min_fdr[ii] <- max(source_fdr[ii], target_fdr[ii])
              }
            }
          }

        }

      }

    }

    background.network$source_score <- source_score
    background.network$target_score <- target_score
    background.network$min_score <- min_score

    background.network$source_fdr <- source_fdr
    background.network$target_fdr <- target_fdr
    background.network$min_fdr <- min_fdr

  }

  return(background.network)

}
