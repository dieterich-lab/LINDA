## ---- message=FALSE-----------------------------------------------------------
library(LINDA)
library(igraph)
library(XML)
library(aggregation)

## ---- warning=FALSE-----------------------------------------------------------
load(file = system.file("extdata", "as_data_toy.RData", package = "LINDA"))
print(as.input)

## ---- warning=FALSE-----------------------------------------------------------
load(file = system.file("extdata", "bg_toy.RData", package = "LINDA"))
print(bg)

## ---- warning=FALSE-----------------------------------------------------------
load(file = system.file("extdata", "tf_act_toy.RData", package = "LINDA"))
print(input.scores)

## ---- message=FALSE-----------------------------------------------------------
res <- runLINDA(input.scores = input.scores, as.input = as.input,
                background.network = bg, solverPath = "~/Downloads/cplex",
                input.node = NULL, pValThresh = 0.05, top = 2, lambda1 = 10,
                lambda2 = 0.001, mipgap = 0.001, relgap = 0.001)

print(res)

## ---- message=FALSE-----------------------------------------------------------
res <- runLINDA(input.scores = input.scores, as.input = NULL,
                background.network = bg, solverPath = "~/Downloads/cplex",
                input.node = NULL, pValThresh = 0.05, top = 2, lambda1 = 10,
                lambda2 = 0.001, mipgap = 0.001, relgap = 0.001)

print(res$combined_interactions)

