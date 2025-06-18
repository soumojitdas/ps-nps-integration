# Bayesian Small Area Estimation: Integrating Probability and Non-Probability Samples

A comprehensive framework for combining probability and non-probability survey samples to improve small area estimation using Bayesian hierarchical models implemented in Stan.

## Project Overview

This repository contains the implementation of a statistical methodology for integrating probability samples (PS) and non-probability samples (NPS) in small area estimation contexts. The project addresses the common challenge where:

- **Probability samples** have known selection probabilities but may observe only proxy measures
- **Non-probability samples** have unknown selection mechanisms but may observe gold standard outcomes
- **Overlap samples** observe both measures when available

Our framework leverages Bayesian hierarchical modeling to produce improved estimates at both population and small area levels.

## Key Features

### Methodological Approaches
- **Hierarchical Bayesian Framework**: Handles small area random effects and uncertainty propagation
- **Multiple Estimation Methods**:
  - OG (Original): Our proposed combined PS-NPS method
  - EMRP: Embedded Multilevel Regression and Poststratification
  - MRP: Standard Multilevel Regression and Poststratification

### Technical Implementation
- **Stan Models**: Efficient Bayesian computation using Hamiltonian Monte Carlo
- **Parallel Processing**: Simulation infrastructure supporting multi-core execution
- **Flexible Sampling Schemes**: SRS, Stratified, and Poisson sampling
- **MAR/MNAR Scenarios**: Handles both Missing at Random and Missing Not at Random mechanisms

## Repository Structure

```
bayesian-small-area-estimation/
├── stan_models/
│   ├── poisson-combined.stan          # Main OG method implementation
│   ├── EMRP-poisson-simp.stan        # EMRP method
│   ├── MRP-N_c.stan                  # MRP method
│   └── stan_model_y2_covariate.stan  # Alternative modeling approach
├── R/
│   ├── simulation/
│   │   ├── empre-re.R                # OG method simulations
│   │   ├── EMPR_predictive.R         # EMRP predictive simulations
│   │   └── MRP_predictive.R          # MRP predictive simulations
│   ├── utils/
│   │   ├── calculate_metrics.R       # Performance metrics
│   │   └── predictive_new.R          # Helper functions
│   └── visualization/
│       └── adv-metrics-plot-mar vs mnar.R  # Results visualization
├── data/
│   └── README.md                     # Data generation instructions
├── results/
│   └── README.md                     # Results storage structure
├── docs/
│   ├── methodology.pdf               # Theoretical framework
│   └── implementation_guide.md       # Usage instructions
└── examples/
    └── simulation_example.R          # Complete workflow example
```

## Getting Started

### Prerequisites
```r
# Required R packages
install.packages(c("cmdstanr", "tidyverse", "future", "progressr"))

# Install cmdstan (for Bayesian computation)
cmdstanr::install_cmdstan()
```

### Basic Usage

```r
# Source required functions
source("R/simulation/empre-re.R")
source("R/utils/calculate_metrics.R")

# Generate population data
pop_data <- generate_population_data(N = 50000)

# Run simulations with OG method
results <- run_OG_simulations(
  S = 100,                    # Number of simulations
  pop_data = pop_data,
  sampling_scheme = "Stratified",
  use_mnar = FALSE,           # MAR scenario
  overlap_pct = 0.11          # 11% overlap
)

# Calculate performance metrics
metrics <- calculate_metrics_unified(
  simulation_results = results,
  true_values = true_values,
  method = "comb"
)
```

## Performance Metrics

The framework evaluates methods using:
- **RMSE**: Root Mean Square Error
- **Bias**: Systematic error in estimates  
- **Coverage**: 95% credible interval coverage
- **Variance Ratio**: Ratio of estimated to true variance

## Simulation Studies

### Sampling Mechanisms
1. **Simple Random Sampling (SRS)**
2. **Stratified Sampling** (by demographic categories)
3. **Poisson Sampling** (probability proportional to size)

### Missing Data Scenarios
- **MAR**: Selection depends only on observed covariates
- **MNAR**: Selection depends on unobserved outcome values

### Small Area Definitions
- Geographic units (states/regions)
- Demographic subgroups
- Cross-classifications

## Key Algorithms

### Algorithm 1: Sample Creation
```
1. Draw probability sample with known weights
2. Draw non-probability sample with bias mechanism  
3. Create controlled overlap
4. Assign observed/missing patterns
```

### Algorithm 2: Estimation Process
```
1. Define cells based on covariates
2. Estimate population cell counts
3. Fit Bayesian hierarchical model
4. Compute estimates
5. Aggregate to small areas
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Implement changes with documentation
4. Add unit tests where applicable
5. Submit a pull request

## Contact

For questions or collaboration opportunities:
- **Email**: soumojit.das5@gmail.com
- **GitHub Issues**: For bug reports and feature requests

---
