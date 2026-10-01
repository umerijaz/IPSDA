# IPSDA: Inverse Particle Swarm Dynamics Analysis for Microbiomes

**IPSDA** is a comprehensive R-based eco-thermodynamic framework designed to decode the non-linear assembly rules of complex microbial ecosystems. By inverting standard Particle Swarm Optimization (PSO), IPSDA maps longitudinal microbiome trajectories onto a continuous kinematic state space. This allows researchers to mathematically extract the physical forces governing community assembly: autecological niche memory ($c_1$) and synecological consensus ($c_2$).

By bridging numerical ecology with statistical mechanics, IPSDA moves beyond correlative networks to calculate ecological escape velocities, forecast catastrophic phase transitions, and map thermodynamic potential energy landscapes in highly dynamic in vivo models.

## Core Features
* **Compositional Data Handling**: Native implementation of robust centered log-ratio (rCLR) transformations to map zero-inflated, constrained simplex data into an unconstrained Euclidean space.
* **Continuous Kinematic Extraction**: Utilizes Generalized Additive Mixed Models (GAMMs) with subject-level random effects to compute continuous velocity and acceleration derivatives from discrete longitudinal sampling.
* **Active-Set Quadratic Programming**: Resolves inverse-PSO equations strictly within non-physical boundaries to prevent mathematically invalid parameter inferences.
* **Thermodynamic Landscaping**: Translates swarm parameters into empirical Kinetic Energy (\(KE\)) and Potential Energy (\(U\)) states to identify systemic vulnerabilities, structural resilience, and dysbiosis tipping points.

## System Requirements & Installation

IPSDA requires R (>= 4.1.0) and relies on several core CRAN and Bioconductor dependencies. 

```R
# Install required CRAN packages
install.packages(c("tidyverse", "mgcv", "nlme", "compositions", "quadprog", "patchwork", "viridis", "ggpubr"))

# Install biomformat from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomformat")
```

## Quick Start & Usage

The core engine is executed via the `run_advanced_ipsda()` wrapper function. The repository includes the `Ijaz2018` longitudinal dataset as a working example. The pipeline natively ingests both `.biom` formats and standard `.csv` abundance matrices.

```R
# Source the IPSDA framework environment
source("IPSDA.R")

# Execute the integrated eco-thermodynamic pipeline on the Ijaz2018 dataset
result <- run_advanced_ipsda(
  count_file         = "Ijaz2018/feature_w_tax.biom",
  meta_file          = "Ijaz2018/meta_table.csv",
  output_dir         = "./ipsda_unified_results",
  sample_id_col      = "SampleID",
  time_col           = "Day", # Continuous time variable
  subject_col        = "Pen_No",  # Controls for repeated measures
  target_vars        = c("BW_Mean", "Gain", "FI", "FCR", "Diet", "Sampling_Groups"), # Metadata mapped to attractors
  window_candidates  = c(3.0, 4.0, 5.0, 6.0, 7.0), # Temporal resolution grid
  window_step        = 1.0,
  inertia_candidates = seq(0.1, 0.8, by = 0.1), # Swarm memory grid
  n_modules          = 4, # Number of synecological guilds
  gamm_k             = 8 # Spline knots for differentiation
)

```

## Key Configuration Parameters

* **`subject_col`**: Essential for robust longitudinal experimental designs. It explicitly informs the GAMM engine of repeated measures, fitting random intercepts to isolate true temporal kinetics from individual host or batch-specific baseline variations.
* **`time_col`**: Must be a continuous numeric vector representing chronological progression (e.g., days, hours, or weeks post-intervention).
* **`target_vars`**: A vector of continuous or categorical metadata covariates. IPSDA automatically regresses the extracted swarm parameters ($c_1$, $c_2$, $KE$, $U$) against these host phenotypes to generate statistical overlays.
* **`window_candidates`** & **`inertia_candidates`**: The framework performs an active-set 2D sliding-window grid search across these temporal spans to maximize local resolution while simultaneously minimizing the variance of the force estimation.
* **`n_modules`**: Dictates the clustering resolution for synecological environmental filtering. Taxa are grouped into co-abundance modules to mathematically model systemic consensus forces.

## Expected Outputs

Upon completion, the specified `output_dir` will populate with raw analytical data structures and high-resolution visualizations:

1. **Optimal Hyperparameter Surfaces**: 2D grid heatmaps evaluating longitudinal model stability and temporal resolution trade-offs.
2. **Phase Space Trajectories**: Visualizations mapping the ecosystem's state shifts between cognitive/autecological ($c_1$) and social/synecological ($c_2$) dominant attractors.
3. **Thermodynamic Energy Landscapes**: Contour maps plotting potential well depth (ecosystem stability) against kinetic volatility (opportunistic colonisation risk).
4. **Escape Velocity Forecasting**: Temporal tracking of the Ecological Escape Index, identifying precise windows where kinetic energy exceeds the local potential well ($KE > \vert{}U\vert{}$), triggering catastrophic phase transitions.
5. **Keystone Driver Extraction**: Ranked physical datasets isolating the specific taxa responsible for establishing deep structural stability versus those driving rapid destabilizing kinetics.

## Citation & Contact

Developed by the **Data Engineering and Systems Biology Lab**, University of Glasgow.

**Corresponding Author**: Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk), Professor of Data Engineering and Systems Biology

**Institution**: James Watt School of Engineering, University of Glasgow

**Repository**: [https://github.com/umerijaz/IPSDA](https://github.com/umerijaz/IPSDA)

If you utilise this computational framework in your research, please cite the corresponding methodology papers covering Inverse Particle Swarm Dynamics Analysis.
