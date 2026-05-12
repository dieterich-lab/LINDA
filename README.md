# LINDA
LINDA R-package - Enio GJERGA

This repository contains the scripts for the ILP (Integer Linear Programming) implementation of the LINDA (Linear Integer programming for Network reconstruction using Transcriptomics and Differential splicing data Analysis) R package and accompanying scripts that implement the method. ILP is a mathematical optimisation algorithm in which the objective function and constraints are linear and the variables are integers. This consensus version supports both hard and soft integration of alternative-splicing (AS) evidence through the `constraints_mode` argument of `runLINDA()`.

### License

Distributed under the GNU GPLv3 License.

### Installation

#### 1. Solver Prerequisites
Before installing LINDA, please keep in mind the following solver pre-requisites:

LINDA requires the interactive version of IBM Cplex as an ILP problem optimiser. The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_|470|135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB).

Once the solvers has been acquired by the user, they must save the executable files in any desired location in the machine they are using and then they can run the LINDA analysis after specifying the solver type (through the *solver* parameter of the functions, 'cplex' or 'cbc'/'lpSolve' in the future) and the path to the executable file (through the *solverPath* parameter).

Soon LINDA is also to run by using the open-source [Coin-Cbc](https://projects.coin-or.org/Cbc) as well as the [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html) R-package for smaller case studies.

**NOTE:** We strongly encourage using CPLEX to solve the LINDA problems since the tool has been mostly maintained by considering CPLEX in mind and also because it showed to be more efficient computationally.

#### 2. Package Dependencies
Additionally before installation, the users must install the following R-package dependencies:
[igraph](https://igraph.org/r/), [aggregation](https://cran.r-project.org/web/packages/aggregation/index.html),
[XML](https://cran.r-project.org/web/packages/XML/index.html) and [aggregation](https://cran.r-project.org/web/packages/aggregation/index.html).

#### 3. Package installation
Once the required solvers have been obtained and the mentioned R-package dependencies have been installed, then the users can proceed with the installation of LINDA.

Currently users can install LINDA directly from the source after downloading the source (tar file) and typing in ```R``` command line the following:

```R
# directly from GitHub
devtools::install_github("dieterich-lab/LINDA", build_vignettes = FALSE)
```

```R
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_LINDA_directory', repos = NULL, type="source")
```

**NOTE:** If you wish for the the Vignettes to be built and for the test example to run successfully, please put the _cplex_ executable file under the /Downloads/ directory and **only** then you can set ```build_vignettes = TRUE```.

## Running LINDA

The LINDA library can be initialized by:

```R
library(LINDA)
```

A documentation of the current version of the main _runLINDA()_ function can be obtained by simply typing ```?runLINDA``` in R once the package has been installed.

## Hard and soft AS integration modes

LINDA can now be run in two AS-aware modes through the `constraints_mode` argument:

- `constraints_mode = "hard"` is the default. AS evidence is converted into explicit ILP constraints. Reactions affected by significant AS events according to `pValThresh` and `splice_effect_sign` are forced to zero, using the same hard constraint formulation as before.
- `constraints_mode = "soft"` keeps those reactions feasible but increases their objective-function penalty according to the AS-derived score. This discourages the solver from selecting AS-affected reactions without making them impossible.

The default `splice_effect_sign` is now `"negative"`, which is appropriate for rMATS-like `IncLevelDifference` values when negative values indicate exon/transcript skipping in the first condition relative to the second condition. You can still set `splice_effect_sign = "positive"` or `splice_effect_sign = "both"`.

Example hard-constrained run:

```R
res_hard <- runLINDA(input.scores = input.scores,
                     as.input = as.input,
                     background.network = bg,
                     solverPath = "~/Downloads/cplex",
                     constraints_mode = "hard",
                     pValThresh = 0.05,
                     splice_effect_sign = "negative",
                     top = 2,
                     lambda1 = 10,
                     lambda2 = 0.001)
```

Example soft-penalty run:

```R
res_soft <- runLINDA(input.scores = input.scores,
                     as.input = as.input,
                     background.network = bg,
                     solverPath = "~/Downloads/cplex",
                     constraints_mode = "soft",
                     splice_effect_sign = "negative",
                     top = 2,
                     lambda1 = 10,
                     lambda2 = 0.001)
```


In a real-case application, depending whether the users are using Transcript Abundance or Exons Skipping data to pinpoint for skipped domains in the analysis, please:

```R
## Load the DIGGER resource with domains mapped to Transcripts:
load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))
```

, or:


```R
## Load the DIGGER resource with domains mapped to Exons:
load(file = system.file("extdata", "digger_human_exons.RData", package = "LINDA"))
```

## Examples

A current example of how to use LINDA over a small Toy test example is documented in the vignettes of the package. For this, please check on the LINDA package [documentation](https://github.com/dieterich-lab/LINDA/blob/main/doc/LINDA.html). Another real case example can be found [here](https://github.com/enio23/LINDA_Example)
