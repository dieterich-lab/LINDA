% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runLINDA.R
\name{runLINDA}
\alias{runLINDA}
\title{Performs the LINDA analysis once the inputs are given.}
\usage{
runLINDA(
  input.scores = input.scores,
  as.input = NULL,
  background.network = bg,
  solverPath = NULL,
  input.node = NULL,
  pValThresh = 0.05,
  pValThresh2 = NULL,
  top = 50,
  lambda1 = 10,
  lambda2 = 0.001,
  mipgap = 0.001,
  relgap = 0.001,
  populate = 1000,
  nSolutions = 100,
  timelimit = 3600,
  intensity = 0,
  replace = 1,
  threads = 0,
  condition = 1,
  solver = "cplex"
)
}
\arguments{
\item{input.scores}{data-frame of two columns containing the estimated
transcription factor activities (TF). The data-frame should have columns ID's
as c('id', 'nes'), indicating to the TF ID's and the corresponding inferred
enrichment score.}

\item{as.input}{a data-frame containing results from the alternative splicing/
exon skipping analysis. It should contain at least three colums with the
following ID's:
id: a transcript/exon ID.
effect: an exon inclusion level difference score / logFC transcript expression
which can be obtained from methods which are used for the detection of
differential AS from replicate RNA-Seq data (i.e. rMATS, MAJIQ, etc.) / methods
for the analysis of differential transcript expression values.
significance: this is as significance value associated with the 'effect'
score. This also can be obtained from differential ES detection tools or from
differential transcript expression analysis.
By default, this parameter is set to NULL which means that LINDA will perform
a splice-un-aware inference of upstream regulatory networks.}

\item{background.network}{data frame of the prior knowledge network of
protein-protein (ppi) and domain-domain (ddi) interactions. The data-frame
should contain at lease 6 columns with the following ID's:
id_source: exon/transcript ID's mapping to the pfam and protein of the source
interactions. In case multiple exons/transcript ID's map to a domain of a
specific protein, then the multiple exons/transcripts are separated by an
underscore symbol - '_'.
id_target: exon/transcript ID's mapping to the pfam and protein of the target
interactions. In case multiple exons/transcript ID's map to a domain of a
specific protein, then the multiple exons/transcripts are separated by an
underscore symbol - '_'.
pfam_source: pfam domain ID's of the interaction sources in the DDI.
pfam_target: pfam domain ID's of the interaction targets in the DDI.
gene_source: gene names of the interaction sources in the PPI.
gene_target: gene names of the interaction targets in the PPI.

In case you using Transcricripts, please load the DIGGER resource with domains
mapped to transcripts:
load(file = system.file("extdata",
                        "digger_human_transcripts.RData",
                        package = "LINDA"))

In case you using Exons, please load the DIGGER resource with domains
mapped to exons:
load(file = system.file("extdata",
                        "digger_human_exons.RData",
                        package = "LINDA"))


}

\item{solverPath}{(optional) location path to the desired solver. By default,
the path to the cplex solver is set to: solverPath = "/usr/bin/cplex"}

\item{input.node}{(optional) a vector of perturbation inputs from where the
protein interactions start. By default input.node=NULL, in which case an
auxilliary 'Perturbation' node is added and which is connected on top of all
the uppermost nodes in the background network. These are the nodes that only
have out-going interactions but no incoming ones.}

\item{pValThresh}{(optional) a p-value threshold to indicate which are the
significantly skipped exons. Those exons in which the corresponding FDR value
(from as.input) falls below the pValThresh parameter, are considered to be
skipped and such information will be integrated in the AS-constraints of the ILP
formulation. By default, pValThresh=0.05.}

\item{splice_effect_sign}{(optional) a character indicating the sign of the
splicing effect to consider for the analysis. It can either be
splice_effect_sign="positive" (for the case when an exon should be considered
as skipped when it crosses the significance threshold pValThresh and it has a
positive sign effect>0); splice_effect_sign="negative" (for the case when
an exon should be considered as skipped when it crosses the significance
threshold pValThresh and it has a negative sign effect<0); or
a positive effect>0 should be considered); splice_effect_sign="negative"
(for the case when only exon skipping event with a negative effect>0
splice_effect_sign="both" (for the case when an exon should be considered as
skipped when it passes the significance threshold pValThresh, regardles of the
sign effect). By default: splice_effect_sign="both".}

\item{top}{(optional) a parameter indicating the number of TF's to be
considered as significantly regulated based on the absolute highest activity
value. The network solution must include as many of the top-regulated TF’s
while penalizing the inclusion of the TF’s which do not appear to be
regulated.By default set to top=50. When considering all TF's, then top="all".}

\item{lambda1}{the penalization term of the primary objective of the objective
function - TF inclusion. This penalty factor is suggested to be set to a
higher value compared to other penalty parameters in order to strongly
penalize the inclusion of not signifcantly regulated TF's. By default,
lambda1=100.}

\item{lambda2}{the penalization terms of the secondary objective of the
objective function - size penalty. The aim of this objective term is to
penalize the inclusion of spurious DDI's in the final solution. By default,
lambda1=1.}

\item{mipgap}{CPLEX parameter which sets an absolute tolerance on the gap
between the best integer objective and the objective of the best node
remaining. When this difference falls below the value of this parameter, the
mixed integer optimization is stopped. By default, mipgap=0.}

\item{relgap}{CPLEX parameter which sets a relative tolerance on the objective
value for the solutions in the solution pool. Solutions that are worse (either
greater in the case of a minimization, or less in the case of a maximization)
than the incumbent solution by this measure are not kept in the solution pool.
For example, if relgap=0.001 (or 0.1 percent), meaning that then solutions worse
than the incumbent by 0.1 percent or more will be discarded. By default,
relgap=0.}

\item{populate}{CPLEX parameter which sets the maximum number of mixed integer
programming (MIP) solutions generated for the solution pool during each call
to the populate procedure. Populate stops when it has generated the amount of
solutions set in this parameter. By default, populate=500.}

\item{nSolutions}{the number of solutions to be provided by LINDA. By default,
nSolutions=100.}

\item{timelimit}{CPLEX parameter which sets the maximum optimization time in
seconds. By default, timelimit=3600.}

\item{intensity}{CPLEX parameter which controls the trade-off between the
number of solutions generated for the solution pool and the amount of time or
memory consumed. Values from 1 to 4 invoke increasing effort to find larger
numbers of solutions. Higher values are more expensive in terms of time and
memory but are likely to yield more solutions. By default, intensity=1.}

\item{replace}{CPLEX parameter which designates the strategy for replacing a
solution in the solution pool when the solution pool has reached its capacity.
The value 0 replaces solutions according to a first-in, first-out policy. The
value 1 keeps the solutions with the best objective values. The value 2
replaces solutions in order to build a set of diverse solutions. By default,
replace=1.}

\item{threads}{CPLEX parameter which manage the number of threads that CPLEX
uses. By default, threads=0 (let CPLEX decide). The number of threads that CPLEX
uses is no more than the number of CPU cores available on the computer where
CPLEX is running.}

\item{condition}{a parameter which can be used in the case when LINDA is
desirde to run over multiple analyses in parallel. It is useful to distinguish
between the multiple ILP problems defined for each case as well as the
solutions obtained. By default, conditions=1.}

\item{solver}{a character parameter indicating the type of solver to be used.
It can only take three different values ('lpSolve', 'cplex' and 'cbc')
depending on the type of solver the user wishes to use. By default,
solver='cplex'.}

\item{solver}{a character parameter indicating the type of solver to be used.
It can only take three different values ('lpSolve', 'cplex' and 'cbc')
depending on the type of solver the user wishes to use. By default,
solver='cplex'.}

}
\value{
Results list containing each of the LINDA unique solutions as well as
the combined ones.
}
\description{
Performs the LINDA analysis once the inputs are given.
}
\examples{
library(LINDA)
library(igraph)
library(XML)

load(file = system.file("extdata", "as_data_toy.RData", package = "LINDA"))
load(file = system.file("extdata", "bg_toy.RData", package = "LINDA"))
load(file = system.file("extdata", "tf_act_toy.RData", package = "LINDA"))

res <- runLINDA(input.scores = input.scores, as.input = as.input, background.network = bg,
solverPath = "~/Downloads/cplex", input.node = NULL, pValThresh = 0.05, top = 2,
lambda1 = 10, lambda2 = 0.001, mipgap = 0.001, relgap = 0.001)

## In case you using Transcricripts, please load the DIGGER resource with domains mapped to transcripts:
# load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))

## In case you using Exons, please load the DIGGER resource with domains mapped to exons:
# load(file = system.file("extdata", "digger_human_exons.RData", package = "LINDA"))

}
