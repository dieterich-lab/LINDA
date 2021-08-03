bin_measurements <- function(input.scores = input.scores, top = 50){

  if(top!="all"){
    if(!is.numeric(top)){
      stop("The top parameter should either be set to 'all' or be numerical")
    } else {
      input.scores$bin <- rep(1, nrow(input.scores))
      input.scores$bin[order(abs(input.scores$nes), decreasing = TRUE)[1:top]] <- -1
    }
  } else {
    input.scores$bin <- rep(-1, nrow(input.scores))
  }

  return(input.scores)

}
