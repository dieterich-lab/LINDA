library(LINDA)
library(igraph)
library(lpSolve)
library(XML)

# get input data
load(file = system.file("extdata", "as_data_toy.RData", package = "LINDA"))
load(file = system.file("extdata", "bg_toy.RData", package = "LINDA"))
load(file = system.file("extdata", "tf_act_toy.RData", package = "LINDA"))

# get expected result
load(file = system.file("result_expected.RData",
                        package="LINDA"))

# obtain actual reesult
result_actual = runLINDA(input.scores = input.scores, as.input = as.input,
                         background.network = bg, input.node = NULL,
                         pValThresh = 0.05, top = 2, solver = "lpSolve")

#testing
test_that("Comparison of the results", {
  expect_equal(result_actual, result_expected)
})
