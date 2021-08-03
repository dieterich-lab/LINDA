integrate_scores_in_bn <- function(as.input = as.input,
                                   background.network = background.network){

  # print("Integrating AS scores in the Background Network..")
  source_score <- rep(0, nrow(background.network))
  target_score <- rep(0, nrow(background.network))
  source_fdr <- rep(1, nrow(background.network))
  target_fdr <- rep(1, nrow(background.network))

  as <- as.input

  uTranscripts <- unique(as$transcript_id)
  for(ii in 1:length(uTranscripts)){

    fdr <- min(as$FDR[which(as$transcript_id==uTranscripts[ii])], na.rm = TRUE)
    score <- min(as$IncLevelDifference[intersect(x = which(as$transcript_id==uTranscripts[ii]),
                                                    y = which(as$FDR==fdr))], na.rm = TRUE)

    idx <- which(grepl(pattern = uTranscripts[ii], x = background.network$exon_source, fixed = TRUE))
    if(length(idx)>0){
      for(jj in 1:length(idx)){
        if(fdr < source_fdr[idx[jj]]){
          source_score[idx[jj]] <- score
          source_fdr[idx[jj]] <- fdr
        }
      }
    }

    idx <- which(grepl(pattern = uTranscripts[ii], x = background.network$exon_target, fixed = TRUE))
    if(length(idx)>0){
      for(jj in 1:length(idx)){
        if(fdr < target_fdr[idx[jj]]){
          target_score[idx[jj]] <- score
          target_fdr[idx[jj]] <- fdr
        }
      }
    }

  }

  min_score <- pmin(source_score, target_score)
  min_fdr <- pmin(source_fdr, target_fdr)

  background.network$source_score <- source_score
  background.network$target_score <- target_score
  background.network$min_score <- min_score

  background.network$source_fdr <- source_fdr
  background.network$target_fdr <- target_fdr
  background.network$min_fdr <- min_fdr

  return(background.network)

}
