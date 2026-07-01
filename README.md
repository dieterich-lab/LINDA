# LINDA

**LINDA** is an R package for **L**inear **IN**teger programming for **D**omain-aware network reconstruction using transcription factor activity and **A**lternative-splicing information.

LINDA reconstructs upstream signalling subnetworks by integrating:

- transcription factor (TF) activity scores,
- a prior-knowledge network of protein-protein and domain-domain interactions,
- optional alternative-splicing (AS) or exon-skipping evidence.

The central function is `runLINDA()`, which supports:

1. **LINDA without AS effects**,
2. **LINDA with hard-constrained AS effects**,
3. **LINDA with soft-constrained AS effects**.

---

## License

Distributed under the GNU GPLv3 License.

---

## Installation

### 1. Solver prerequisites

LINDA formulates the reconstruction problem as an integer linear programming (ILP) problem and therefore requires an ILP solver.

The recommended solver is **IBM ILOG CPLEX**, especially for larger analyses. After installation, provide the path to the CPLEX executable through `solverPath`, for example:

```r
solverPath = "~/Downloads/cplex"
```

### 2. Install R dependencies

```r
install.packages(c(
  "devtools",
  "igraph",
  "XML",
  "aggregation",
  "readxl",
  "readr",
  "rmarkdown",
  "qpdf"
))
```

### 3. Install LINDA

From GitHub (development branch):

```r
devtools::install_github("https://github.com/dieterich-lab/LINDA", build_vignettes = FALSE)
```

From local source:

```r
install.packages("path_to_extracted_LINDA_directory", repos = NULL, type = "source")
```

---

## Load LINDA

```r
library(LINDA)
library(igraph)
library(XML)
library(aggregation)
```

Open the help page with:

```r
?runLINDA
```

---

## Main function

```r
runLINDA(
  input.scores,
  as.input = NULL,
  background.network,
  solverPath = NULL,
  input.node = NULL,
  constraints_mode = "hard",
  pValThresh = 0.05,
  splice_effect_sign = "negative",
  top = 50,
  lambda1 = 100,
  lambda2 = 1,
  mipgap = 0,
  relgap = 0,
  solver = "cplex"
)
```

### Core inputs

#### `input.scores`
A data frame with TF activity values. It should contain at least:

- `id`: TF identifier,
- `nes`: activity or enrichment score.

#### `as.input`
A data frame with alternative-splicing information. It should contain at least:

- `id`: transcript or exon identifier,
- `effect`: splicing effect size,
- `significance`: associated significance value.

If `as.input = NULL`, LINDA runs **without AS effects**.

#### `background.network`
A background network describing domain-aware directed interactions, typically with fields such as:

- `exon_source`, `exon_target`,
- `pfam_source`, `pfam_target`,
- `gene_source`, `gene_target`.

---

## AS integration modes

### 1. LINDA without AS effects
Use:

```r
as.input = NULL
```

In this case, no splice information is used, and `constraints_mode` has no practical effect on the result.

### 2. LINDA with hard-constrained AS effects
Use:

```r
constraints_mode = "hard"
```

This is the default mode. Significant AS-affected interactions are treated as strict constraints and are excluded from the feasible solution.

### 3. LINDA with soft-constrained AS effects
Use:

```r
constraints_mode = "soft"
```

In this mode, AS-affected interactions are penalized rather than strictly removed, so the solver may still use them if needed.

### Direction of AS effect

The `splice_effect_sign` argument controls which AS direction is treated as functionally relevant:

```r
splice_effect_sign = "negative"
splice_effect_sign = "positive"
splice_effect_sign = "both"
```

The default is:

```r
splice_effect_sign = "negative"
```

---

# Toy example

This section documents the built-in toy example step by step and organizes it into the three analysis modes:

1. **LINDA without AS effects**,
2. **LINDA with hard-constrained AS effects**,
3. **LINDA with soft-constrained AS effects**.

The example uses:

- a toy background network,
- a toy AS-effect table,
- toy TF activity scores.

## Step 1: Load required packages

```r
library(LINDA)
library(igraph)
library(XML)
library(aggregation)
```

## Step 2: Define the solver path

```r
solver_path <- "~/Downloads/cplex"
```

Replace this path with the location of the CPLEX executable on your machine.

## Step 3: Load the toy input objects

```r
load(file = system.file("extdata", "as_data_toy.RData", package = "LINDA"))
print(as.input)

load(file = system.file("extdata", "bg_toy.RData", package = "LINDA"))
print(bg)

load(file = system.file("extdata", "tf_act_toy.RData", package = "LINDA"))
print(input.scores)
```

## Toy-example overview

The figure below summarizes the toy data used in the three example analyses. It shows:

- the background network,
- the AS-effect table,
- the TF activity values,
- the toy domains/transcripts highlighted as splice-affected in this example.

![Toy-example inputs and overview](man/figures/toy_example_info.jpg)

**Interpretation of the toy setup**

- **P1–P5** denote proteins/nodes in the signalling network.
- **D1–D14** denote domains or domain-associated transcript/exon features.
- **M1–M3** denote transcription factor readouts.
- The AS-effect table marks which transcript/domain-associated elements are significantly splice-affected.
- In this toy example, `top = 2` means that LINDA focuses on the two strongest TF signals according to absolute activity value.

---

# i) LINDA without AS effects

When no AS input is provided, LINDA reconstructs the signalling network using only the background network and TF activity. This is the same regardless of whether `constraints_mode` is set to `"hard"` or `"soft"`, because no AS information is available to constrain or penalize the solution.

## Code

```r
res_no_as <- runLINDA(
  input.scores = input.scores,
  as.input = NULL,
  background.network = bg,
  solverPath = solver_path,
  input.node = NULL,
  top = 2,
  lambda1 = 10,
  lambda2 = 0.001,
  mipgap = 0,
  relgap = 0
)

print(res_no_as$combined_interactions)
```

## What this run does

- ignores alternative-splicing information,
- uses the top 2 TF activity readouts,
- identifies the shortest / most parsimonious signalling solution that explains the selected TFs.

## Visualization

The figure below illustrates the splice-unaware solution. Since no AS effect is considered, LINDA selects the shortest signalling paths needed to reach the relevant TF outputs.

![LINDA without AS effects](man/figures/toy_example_no_as.jpg)

**Key point:** This is the baseline reference solution against which the AS-aware hard and soft modes can be compared.

---

# ii) LINDA with hard-constrained effects

In this mode, significant splice-affected interactions are treated as hard constraints. With the toy example shown here, splice-affected negative transcript regulation together with a significance threshold of `0.05` causes certain connections to be skipped.

## Code

```r
res_hard <- runLINDA(
  input.scores = input.scores,
  as.input = as.input,
  background.network = bg,
  solverPath = solver_path,
  constraints_mode = "hard",
  pValThresh = 0.05,
  input.node = NULL,
  top = 2,
  lambda1 = 10,
  lambda2 = 0.001,
  mipgap = 0,
  relgap = 0
)

print(res_hard$combined_interactions)
```

## What this run does

- includes the AS input,
- applies AS information as strict exclusion constraints,
- removes splice-affected edges from the feasible solution space.

## Visualization

The figure below illustrates the hard-constrained solution.

![LINDA with hard-constrained AS effects](man/figures/toy_example_hard_as.jpg)

**Interpretation:** Connections involving the splice-affected elements (illustrated in the toy example as D3_P1, D7_P3 and D13_P5 in the figure annotation) are skipped, so LINDA is forced to identify an alternative feasible route.

---

# iii) LINDA with soft-constrained effects

In this mode, splice-affected interactions are not absolutely forbidden. Instead, they are penalized in the objective function, so LINDA prefers cleaner alternatives when available, but may still use splice-affected edges if they are required to explain the TF pattern.

## Code

```r
res_soft <- runLINDA(
  input.scores = input.scores,
  as.input = as.input,
  background.network = bg,
  solverPath = solver_path,
  constraints_mode = "soft",
  input.node = NULL,
  top = 2,
  lambda1 = 10,
  lambda2 = 0.001,
  mipgap = 0,
  relgap = 0
)

print(res_soft$combined_interactions)
```

## What this run does

- includes the AS input,
- penalizes splice-affected edges instead of removing them,
- allows LINDA to trade off path optimality against splice penalties.

## Visualization

The figure below illustrates the soft-constrained solution.

![LINDA with soft-constrained AS effects](man/figures/toy_example_soft_as.jpg)

**Interpretation:** Connections involving splice-affected elements are discouraged rather than blocked. Therefore, LINDA may prefer a route such as **P2 → P4** over **P1 → P3** if the overall splice penalty is lower.

---

## Complete toy-example script

```r
library(LINDA)
library(igraph)
library(XML)
library(aggregation)

solver_path <- "~/Downloads/cplex"

load(file = system.file("extdata", "as_data_toy.RData", package = "LINDA"))
load(file = system.file("extdata", "bg_toy.RData", package = "LINDA"))
load(file = system.file("extdata", "tf_act_toy.RData", package = "LINDA"))

# i) LINDA without AS effects
res_no_as <- runLINDA(
  input.scores = input.scores,
  as.input = NULL,
  background.network = bg,
  solverPath = solver_path,
  input.node = NULL,
  top = 2,
  lambda1 = 10,
  lambda2 = 0.001,
  mipgap = 0,
  relgap = 0
)
print(res_no_as$combined_interactions)

# ii) LINDA with hard-constrained effects
res_hard <- runLINDA(
  input.scores = input.scores,
  as.input = as.input,
  background.network = bg,
  solverPath = solver_path,
  constraints_mode = "hard",
  pValThresh = 0.05,
  input.node = NULL,
  top = 2,
  lambda1 = 10,
  lambda2 = 0.001,
  mipgap = 0,
  relgap = 0
)
print(res_hard$combined_interactions)

# iii) LINDA with soft-constrained effects
res_soft <- runLINDA(
  input.scores = input.scores,
  as.input = as.input,
  background.network = bg,
  solverPath = solver_path,
  constraints_mode = "soft",
  input.node = NULL,
  top = 2,
  lambda1 = 10,
  lambda2 = 0.001,
  mipgap = 0,
  relgap = 0
)
print(res_soft$combined_interactions)
```

---

## Understanding the output

Each `runLINDA()` call returns a list object. The most useful elements are:

```r
res$separate_interactions
res$combined_interactions
```

### `res$separate_interactions`
Contains the individual LINDA solutions found by the solver.

### `res$combined_interactions`
Contains the merged or summarized interaction set across the returned solutions and is usually the most convenient object for downstream inspection.

For example:

```r
print(res_hard$combined_interactions)
View(res_hard$combined_interactions)
```

---

## Suggested comparison of the three modes

```r
res_no_as$combined_interactions
res_hard$combined_interactions
res_soft$combined_interactions
```

A useful interpretation workflow is:

- compare **without AS** vs **hard-constrained** to see which edges are fully removed by splice evidence,
- compare **without AS** vs **soft-constrained** to see how splice penalties redirect the chosen path,
- compare **hard-constrained** vs **soft-constrained** to assess the difference between strict exclusion and penalized preference.

---

## Real-case DIGGER resources

For real-case analyses, LINDA requires a background network in which protein and domain interactions are mapped to transcript-level or exon-level identifiers. The package includes both original DIGGER-derived resources and DIGGER v2-derived resources.

The original DIGGER resource was developed to explore how alternative splicing can affect protein-protein interactions by integrating protein-protein interactions, domain-domain interactions and residue-level interaction information. DIGGER 2.0 extends this idea and adds support for both human and mouse resources, including experimentally supported and predicted domain-domain interactions.

### Original DIGGER resources

The original [DIGGER](https://academic.oup.com/nar/article/49/D1/D309/5911747) resources can be loaded as follows:

```r
load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))
load(file = system.file("extdata", "digger_human_exons.RData", package = "LINDA"))
load(file = system.file("extdata", "digger_human.RData", package = "LINDA"))
load(file = system.file("extdata", "digger_mouse.RData", package = "LINDA"))
```

### DIGGER v2 resources

LINDA also includes [DIGGER v2-based](https://academic.oup.com/nar/article/53/W1/W245/8126897) background networks for both **human** and **mouse**. These resources are provided at two mapping resolutions:

- **exon-level** resources,
- **transcript-level** resources.

For each organism and resolution, the DIGGER v2 resources are split into three confidence classes:

- **Gold**: highest-confidence interaction class,
- **Silver**: intermediate-confidence interaction class,
- **Bronze**: broader/lower-confidence interaction class.

The available DIGGER v2 resource files are:

| Organism | Resolution | Confidence | File |
|---|---|---:|---|
| Human | Exons | Gold | `digger_v2_human_exons_gold.RData` |
| Human | Exons | Silver | `digger_v2_human_exons_silver.RData` |
| Human | Exons | Bronze | `digger_v2_human_exons_bronze.RData` |
| Human | Transcripts | Gold | `digger_v2_human_transcripts_gold.RData` |
| Human | Transcripts | Silver | `digger_v2_human_transcripts_silver.RData` |
| Human | Transcripts | Bronze | `digger_v2_human_transcripts_bronze.RData` |
| Mouse | Exons | Gold | `digger_v2_mouse_exons_gold.RData` |
| Mouse | Exons | Silver | `digger_v2_mouse_exons_silver.RData` |
| Mouse | Exons | Bronze | `digger_v2_mouse_exons_bronze.RData` |
| Mouse | Transcripts | Gold | `digger_v2_mouse_transcripts_gold.RData` |
| Mouse | Transcripts | Silver | `digger_v2_mouse_transcripts_silver.RData` |
| Mouse | Transcripts | Bronze | `digger_v2_mouse_transcripts_bronze.RData` |

### Loading one DIGGER v2 resource

For example, to load the human exon-level gold-confidence DIGGER v2 network:

```r
load(file = system.file(
  "extdata",
  "digger_v2_human_exons_gold.RData",
  package = "LINDA"
))
```

If you want to load resources programmatically without depending on the internal object name stored inside the `.RData` file, you can use a small helper function:

```r
load_linda_extdata <- function(filename) {
  env <- new.env(parent = emptyenv())
  load(system.file("extdata", filename, package = "LINDA"), envir = env)
  object_names <- ls(env)

  if (length(object_names) != 1) {
    stop("Expected exactly one object in ", filename,
         "; found: ", paste(object_names, collapse = ", "))
  }

  env[[object_names[1]]]
}

bg_human_exons_gold <- load_linda_extdata("digger_v2_human_exons_gold.RData")
```

### Combining confidence classes

If you want to combine different confidence classes, load the desired resources and combine them with `rbind()`.

For example, to combine the human exon-level gold, silver and bronze networks:

```r
bg_human_exons_gold <- load_linda_extdata("digger_v2_human_exons_gold.RData")
bg_human_exons_silver <- load_linda_extdata("digger_v2_human_exons_silver.RData")
bg_human_exons_bronze <- load_linda_extdata("digger_v2_human_exons_bronze.RData")

bg_human_exons_all <- rbind(
  bg_human_exons_gold,
  bg_human_exons_silver,
  bg_human_exons_bronze
)
```

If you want to keep track of the confidence class after combining, add a class label before using `rbind()`:

```r
bg_human_exons_gold$confidence_class <- "gold"
bg_human_exons_silver$confidence_class <- "silver"
bg_human_exons_bronze$confidence_class <- "bronze"

bg_human_exons_all <- rbind(
  bg_human_exons_gold,
  bg_human_exons_silver,
  bg_human_exons_bronze
)
```

You can use the same approach for mouse resources and for transcript-level resources. For example:

```r
bg_mouse_transcripts_gold <- load_linda_extdata("digger_v2_mouse_transcripts_gold.RData")
bg_mouse_transcripts_silver <- load_linda_extdata("digger_v2_mouse_transcripts_silver.RData")

bg_mouse_transcripts_gold_silver <- rbind(
  bg_mouse_transcripts_gold,
  bg_mouse_transcripts_silver
)
```

Real-case applications of LINDA over the [ENCORE Initiative](https://www.encodeproject.org/encore-matrix/?type=Experiment&status=released&internal_tags=ENCORE) dataset by using the DIGGER-v1 background have been provided in this separate [GitHub repo](https://github.com/enio23/LINDA_ENCORE).

### References

Gjerga et al., **Characterizing alternative splicing effects on protein interaction networks with LINDA**, *Bioinformatics*, 2023. <https://doi.org/10.1093/bioinformatics/btad224>
