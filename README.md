# LINDA
LINDA R-package - Enio GJERGA

This repository contains the scripts for the ILP (Integer Linear Programming) implementation of the LINDA R package and accompanying scripts that implement the method. ILP is a mathematical optimisation algorithm in which the objective function and constraints are linear and the variables are integers.

### License

Distributed under the GNU GPLv3 License.

### Installation

#### 1. Solver Prerequisites
Before installing LINDA, please keep in mind the following solver pre-requisites:

LINDA requires the interactive version of IBM CPLEX as an ILP problem optimiser. The IBM ILOG CPLEX is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_|470|135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB). With such a solver, users are bale to train large-scale models and enumerate multiple alternative solutions.

Alternatively to CPLEX, users can rely on the COIN-OR branch and cut (or Cbc) solver. This solver is an **open-source mixed integer programming solver** written in C++ and is freely available [here](https://projects.coin-or.org/Cbc). Similar to CPLEX, by sing the CBC solver we are able to train large-scale models however users are not able to enumerate more than 1 solution.

Additionally, users can rely on the free open-source [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html) solver which has been implemented as a CRAN package. Such a solver however is only suitable to be used for small Toy case-studies.

Once the solvers has been acquired by the user, they must save the executable files in any desired location in the machine they are using and then they can run the LINDA analysis after specifying the solver type (through the *solver* parameter of the functions 'cplex', 'cbc' or 'lpSolve') and the path to the executable file (through the *solverPath* parameter for the 'cplex' or 'cbc' case).

**NOTE:** We strongly encourage using cplex to solve the LINDA problems since the tool has been mostly maintained by considering CPLEX in mind and also because it showed to be more efficient computationally.

#### 2. Package Depedencies
Additionally before installation, the users must install the following R-package depedencies:
[igraph](https://igraph.org/r/) and
[XML](https://cran.r-project.org/web/packages/XML/index.html).

#### 3. Package installation
Once the required solvers have been obtained and the mentioned R-package depedencies have been installed, then the users can proceed with the installation of LINDA.

Currently users can install LINDA directly from the source after downloading the source (tar file) and typing in ```R``` command line the following:

```R
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_LINDA_directory', repos = NULL, type="source")
```

Otherwise, it is possible to install directly form github:
```R
install.packages("devtools")
devtools::install_github("dieterich-lab/LINDA")
```

## Running LINDA

The LINDA library can be initialized by:

```R
library(LINDA)
```

A documentation of the current version of the main _runLINDA()_ function can be obtained by simply typing ```?runLINDA``` in R once the package has been installed.

## Examples

A current example of how to use LINDA over a small Toy test example is ducumented in the vignettes of the package. For this, please check on the LINDA package [documentation](https://github.com/dieterich-lab/LINDA/blob/main/doc/LINDA.html).

