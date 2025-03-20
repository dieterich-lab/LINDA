integrate_scores_in_bn_soft <- function(as.input = as.input,
                                        background.network = background.network,
                                        splice_effect_sign = splice_effect_sign){
  
  if(!is.null(as.input)){
    colnames(as.input) <- c("exon_id", "IncLevelDifference", "pval")
  }
  
  idx <- which(as.input$pval==0)
  if(length(idx)>0){
    as.input$pval[idx] <- 0.000001
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
    
    for(ii in seq_len(nrow(background.network))){
      
      # source
      transcripts <- unique(unlist(strsplit(x = background.network$exon_source[ii], split = "_", fixed = TRUE)))
      idx <- which(as$exon_id%in%transcripts)
      if(length(idx)>0){
        
        source_score[ii] <- mean(as$IncLevelDifference[idx])
        source_fdr[ii] <- fisher(as$pval[idx])
        
      }
      
      # target
      transcripts <- unique(unlist(strsplit(x = background.network$exon_target[ii], split = "_", fixed = TRUE)))
      idx <- which(as$exon_id%in%transcripts)
      if(length(idx)>0){
        
        target_score[ii] <- mean(as$IncLevelDifference[idx])
        target_fdr[ii] <- fisher(as$pval[idx])
        
      }
      
    }
    
    if(splice_effect_sign == "negative"){
      
      min_score <- rep(0, nrow(background.network))
      min_fdr <- rep(1, nrow(background.network))
      for(ii in 1:nrow(background.network)){
        
        if(source_score[ii]<0 || target_score[ii]<0){
          
          min_score[ii] <- (-1)*log2(min(c(source_fdr[ii], target_fdr[ii])))
          min_fdr[ii] <- min(c(source_fdr[ii], target_fdr[ii]))
          
        }
        
      }
      
    } else {
      
      if(splice_effect_sign == "positive"){
        
        min_score <- rep(0, nrow(background.network))
        min_fdr <- rep(1, nrow(background.network))
        for(ii in 1:nrow(background.network)){
          
          if(source_score[ii]>0 || target_score[ii]>0){
            
            min_score[ii] <- (-1)*log2(min(c(source_fdr[ii], target_fdr[ii])))
            min_fdr[ii] <- min(c(source_fdr[ii], target_fdr[ii]))
            
          }
          
        }
        
      } else {
        
        if(splice_effect_sign == "both"){
          
          min_score <- rep(0, nrow(background.network))
          min_fdr <- rep(1, nrow(background.network))
          
          for(ii in 1:nrow(background.network)){
            
            min_score[ii] <- (-1)*log2(min(c(source_fdr[ii], target_fdr[ii])))
            min_fdr[ii] <- min(c(source_fdr[ii], target_fdr[ii]))
            
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
