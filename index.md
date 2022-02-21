## Overview

LINDA (Linear Integer programming for Network reconstruction using transcriptomics and Differential splicing data Analysis), is an R-Package used to identify the mechanistic upstream regulatory signalling processes that drive gene expression changes by also taking into account the contribution of alternative splicing events to signal protein variability and information flow modulation.

### Install

LINDA can be installed either locally by downloading its [source code](https://github.com/dieterich-lab/LINDA) or directly from GitHub via [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html).

#### Download
1.  Download LINDA package (*Main* or *Development* branch) by clicking on **Code** and then **Download ZIP**.
2.  Start R.
3.  Unzip package and set working directory to where LINDA has been downloaded in the R workspace.
4.  Install LINDA by running: `install.packages('LINDA/', repos = NULL, type="source")`.

#### devtools
1.  Make sure [devtools](https://github.com/dieterich-lab/LINDA) is installed.
2.  Start R.
3.  Run `install_github("dieterich-lab/LINDA", build_vignettes = TRUE)`.

### Documentation
The LINDA pipeline and it prerequisites have been documented on its Vignettes. Also documentation about to test cases (one Toy example and one Real-case application) have been provided.
1.  For an explanation of the LINDA pipeline and a to run a small Toy test case-study, you can check the vignettes by simply running: `vignette("LINDA")`.
2.  Also a step-by-step explanation of LINDA application over a real case-study has been provided in: [https://github.com/enio23/LINDA_Example](https://github.com/enio23/LINDA_Example).

### Support or Contact

Having trouble with LINDA? You can open an Issue on the dedicated [LINDA repository](https://github.com/dieterich-lab/LINDA), or simply drop us an email on [E.Gjerga@uni-heidelberg.de](E.Gjerga@uni-heidelberg.de) & [Christoph.Dieterich@uni-heidelberg.de](Christoph.Dieterich@uni-heidelberg.de) and we’ll help you sort it out.
