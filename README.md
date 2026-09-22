# Length-Structured Population Model for Peruvian Anchoveta

This repository contains an R implementation of a **length-structured population projection model** developed for the north-central stock of Peruvian anchoveta (*Engraulis ringens*).

The model projects abundance at length on a weekly time step, incorporating growth, natural mortality, catch removals, maturity, and biomass. It also estimates fishing mortality by length class and week using the Baranov catch equation and derives a set of fishing mortality indicators.

## Repository structure

```text
.
├── code/
│   ├── 0_internal_function.R
│   ├── 1_model_functions.R
│   └── 2_run_model.R
│
└── input/
    ├── surveyFile.RData
    ├── catchFile.RData
    ├── life_history_params_by_envir_scenario.RData
    └── species_data.RData
```

### `0_internal_function.R`

Contains the internal functions used by the population model, including:

* growth calculations;
* length-transition matrices;
* natural mortality at length;
* population projection;
* natural mortality by length class;
* fishing mortality estimation using the Baranov equation;
* fishing mortality summaries and indicators.

### `1_model_functions.R`

Contains the main functions used to run the model:

* `projectmodel()` — projects abundance, biomass, and mature biomass for one or multiple abundance replicates;
* `lengthBasedModel()` — performs the weekly population projection across bootstrap replicates and summarizes abundance and biomass.

### `2_run_model.R`

Main script used to:

1. load input data;
2. define model settings;
3. prepare catch-at-length and survey abundance-at-length matrices;
4. run the weekly population projection;
5. estimate fishing mortality;
6. calculate fishing mortality indicators;
7. store the model outputs.

---

## Requirements

The model requires R and the following packages:

```r
install.packages(c("dplyr", "progress"))
```

The scripts use:

```r
library(dplyr)
```

The `progress` package is called internally using `progress::progress_bar()`.

---

## Input data

Four `.RData` files are required.

### Survey abundance

```text
input/surveyFile.RData
```

The file should contain the object:

```r
survey
```

representing survey abundance by length class. Columns can contain bootstrap replicates of abundance-at-length.

Survey abundance is assumed to be expressed in **millions of individuals** in the current model configuration.

---

### Catch data

```text
input/catchFile.RData
```

The file should contain:

```r
catchFile
```

with weekly catch-at-length information.

The current implementation expects length-class columns corresponding to:

```text
2.0, 2.5, 3.0, ..., 19.5, 20.0 cm
```

and uses observations where:

```r
type == "Week"
```

to construct the weekly catch-at-length matrix.

Catch is expressed as **number of individuals** in the current configuration.

---

### Environmental scenarios

```text
input/life_history_params_by_envir_scenario.RData
```

The file should contain the object:

```r
scenarios
```

Each environmental scenario contains the biological parameters used by the model, including:

* `k` — growth coefficient;
* `Linf` — asymptotic length;
* `sizeM` — length thresholds used to define natural mortality groups;
* `vectorM` — natural mortality values associated with those groups.

The environmental scenario used in the projection is selected through:

```r
environmental_scenario <- "neutral"
```

---

### Species parameters

```text
input/species_data.RData
```

The file should contain:

```r
species_data
```

with the biological information required to define the length structure and maturity ogive.

The internal functions use the following parameters:

```text
Lmin
Lmax
bin
mat1
mat2
```

The length-weight relationship is defined as:

```text
W = a × L^b
```

with the current values:

```r
a = 0.0072
b = 2.9855
```

---

## Running the model

Set the repository root as the R working directory and run:

```r
source("code/2_run_model.R")
```

The main script automatically loads the required functions:

```r
source("code/0_internal_function.R")
source("code/1_model_functions.R")
```

and the input `.RData` files.

---

## Model workflow

The model follows the sequence:

```text
Survey abundance-at-length
          │
          ▼
Initial abundance
          │
          ▼
Natural mortality
(first half of week)
          │
          ▼
Catch removal
(midpoint of week)
          │
          ▼
Natural mortality
(second half of week)
          │
          ▼
Growth / length transition
          │
          ▼
Abundance-at-length
for the next week
```

This process is repeated for every fishing week and for every bootstrap replicate of the initial survey abundance.

Growth is represented through a length-transition matrix based on the relationship:

```text
L(t + Δt) = Linf - [Linf - L(t)] exp(-k Δt)
```

where growth parameters depend on the selected environmental scenario.

Catch is removed at the midpoint of each weekly interval. Natural mortality is therefore applied in two half-steps, before and after catch removal.

When projected catch exceeds the available abundance in a length class, negative abundance values are truncated to zero and a warning is generated.

---

## Bootstrap projection

When:

```r
data_type <- "boot"
```

the model propagates all survey bootstrap replicates independently.

For each time step, the median abundance across bootstrap replicates is calculated.

The returned object contains both the median trajectories and the complete bootstrap arrays.

---

## Fishing mortality

Fishing mortality \(F\) is estimated independently for each length class and fishing week using the Baranov catch equation.

The numerical solution is obtained using `uniroot()`.

The estimation requires abundance and catch to be expressed in compatible units.

The model subsequently calculates several fishing mortality indicators.

### Adult fishing mortality

Abundance-weighted fishing mortality for individuals with:

```r
length >= adult_length
```

The default threshold is:

```r
adult_length = 12
```

### Juvenile fishing mortality

Abundance-weighted fishing mortality for individuals below the adult length threshold.

### Stock fishing mortality

Abundance-weighted fishing mortality calculated across all length classes.

### Maximum weekly fishing mortality

Maximum estimated \(F\) among length classes for each week.

### Observed selectivity

Relative fishing mortality-at-length, normalized by the maximum value:

```text
selectivity = mean F-at-length / maximum mean F-at-length
```

### Cumulative fishing mortality

Fishing mortality accumulated across fishing weeks for each length class.

---

## Model output

The main script stores the results in:

```r
results_model
```

with the following structure:

```r
results_model <- list(
  CatchData   = CatchData,
  SurveyData  = SurveyData,
  model       = output,
  M_at_length = M_l,
  F_at_length = F_model,
  indicators  = indicadores
)
```

### `results_model$model`

Contains the population projection, including:

```text
modelN
```

Median abundance-at-length through time.

```text
modelB
```

Median biomass-at-length through time.

```text
Biomass
```

Total biomass through time.

When bootstrap output is requested, the object also contains:

```text
modelArrayN
modelArrayB
```

with the complete bootstrap distributions.

### `results_model$M_at_length`

Natural mortality associated with each length class.

### `results_model$F_at_length`

Fishing mortality matrix where rows represent length classes and columns represent fishing weeks.

### `results_model$indicators`

Contains:

```text
F_adults
F_juveniles
F_stock
F_max
selectivity
cumulative_F
```

---

## Current model configuration

The example implementation is configured for:

```text
Stock: North-Central Peruvian anchoveta
Species: Engraulis ringens
Length range: 2–20 cm
Length interval: 0.5 cm
Time step: weekly
Temporal frequency: 52 steps per year
Adult length threshold: 12 cm
Default environmental scenario: neutral
```

Survey abundance is expressed in millions of individuals and internally converted to number of individuals before population projection.

---

## Notes

The model assumes that:

* population dynamics are represented by discrete length classes;
* growth can be represented by a deterministic length-transition matrix;
* biological parameters can vary among predefined environmental scenarios;
* natural mortality varies among length groups;
* catch removal occurs at the midpoint of the weekly interval;
* survey abundance provides the initial abundance-at-length distribution;
* uncertainty in initial abundance can be propagated through bootstrap replicates.

---

## Example outputs

After running the model, total projected biomass can be inspected with:

```r
results_model$model$Biomass
```

and weekly stock fishing mortality with:

```r
results_model$indicators$F_stock
```
