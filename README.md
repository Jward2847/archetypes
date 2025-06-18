# The epidemiology of pathogens with pandemic potential: A review of key parameters and clustering analysis 

This project implements a Monte Carlo approach combined with K-means clustering to identify distinct pathogen archetypes based on epidemiological parameters. The analysis accounts for parameter uncertainty through probabilistic sampling and provides robust clustering results through ensemble methods.

## Project Overview

The goal is to classify infectious disease pathogens into distinct archetypes based on their epidemiological characteristics, including:
- **R0 (Basic Reproduction Number)**: Transmission potential
- **Serial Interval (SI)**: Time between symptom onset in successive cases
- **Case Fatality Rate (CFR)**: Proportion of cases that result in death
- **Presymptomatic Transmission Proportion**: Fraction of transmission occurring before symptom onset
- **Transmission Route**:
  
## Scripts Overview

### 1. `mcmc_Kmeans_analysis.R` - Main Analysis Script

**Purpose**: Performs the complete MCMC K-means clustering analysis with K=6 clusters.

**Key Features**:
- **MCMC Sampling**: Generates 5,000 parameter samples per pathogen using uncertainty distributions
- **Per-iteration Clustering**: Applies K-means clustering to each MCMC iteration
- **Ensemble/Consensus Clustering**: Combines results across all iterations using hierarchical clustering
- **Comprehensive Output**: Produces cluster assignments, centroids, and visualisations

**Main Steps**:
1. **Parameter Sampling**: Samples from uncertainty distributions (normal, beta, gamma, lognormal, uniform)
2. **Presymptomatic Proportion Estimation**: Calculates from serial interval and incubation period distributions
3. **Per-iteration K-means**: Clusters pathogens for each MCMC iteration (K=6)
4. **Co-assignment Matrix**: Tracks how often pathogens are grouped together
5. **Consensus Clustering**: Uses hierarchical clustering on dissimilarity matrix
6. **Visualisation**: PCA plots and dendrograms with full pathogen names

**Outputs**:
- `mcmc_parameter_samples.csv`: All MCMC parameter samples
- `consensus_cluster_assignments_k6.csv`: Final cluster assignments
- `cluster_centroids_summary_with_ci_k6.csv`: Cluster characteristics with confidence intervals
- `consensus_dendrogram_colored_k6.png`: Hierarchical clustering dendrogram
- `pca_modal_cluster_plot_k6.png`: PCA visualisation of clusters

### 2. `find_optimal_K.R` - Optimal Cluster Number Analysis

**Purpose**: Determines the optimal number of clusters (K) using the silhouette method.

**Key Features**:
- **Silhouette Analysis**: Tests K values from 2 to 10
- **Dissimilarity-based Evaluation**: Uses the consensus clustering dissimilarity matrix
- **Visual Results**: Plots average silhouette width vs. K
- **Optimal K Selection**: Identifies K with maximum silhouette width

**Process**:
1. Loads the dissimilarity matrix from the main analysis
2. Performs hierarchical clustering with average linkage
3. Cuts the tree at different K values (2-10)
4. Calculates silhouette width for each K
5. Identifies optimal K based on maximum silhouette width

**Outputs**:
- `optimal_k_silhouette_plot.png`: Plot showing silhouette width vs. K
- `optimal_k_silhouette_results.csv`: Numerical results for each K value

### 3. `mcmc_sensitivity_analyis.R` - Sensitivity Analysis

**Purpose**: Tests the robustness of clustering results by varying included pathogens and parameters.

**Two Sensitivity Analyses**:

#### RVFV Sensitivity Analysis
- **Data**: Includes RVFV pathogen (`pathogen_params_RVFV.csv`)
- **Excluded Parameters**: Serial interval and presymptomatic proportion
- **Clustering Variables**: R0, CFR, and transmission routes
- **Optimal K**: 6 (determined by silhouette analysis)

#### HIV Sensitivity Analysis  
- **Data**: Includes HIV pathogen (`pathogen_params_HIV.csv`)
- **Excluded Parameters**: Serial interval and case fatality rate
- **Clustering Variables**: R0, presymptomatic proportion, and transmission routes
- **Optimal K**: 8 (determined by silhouette analysis)

**Key Features**:
- **Conditional Parameter Sampling**: Only samples parameters needed for each analysis
- **Separate Output Directories**: Results saved to `MCMC_sens_analysis/`
- **Comprehensive Evaluation**: Includes silhouette analysis for each scenario
- **Robustness Assessment**: Compares clustering stability across different parameter sets

## Data Requirements

### Input Files
- `pathogen_params_kmeans.csv`: Main parameter file with 19 pathogens
- `pathogen_params_RVFV.csv`: Extended dataset including RVFV
- `pathogen_params_HIV.csv`: Extended dataset including HIV

### Parameter Structure
Each pathogen requires:
- **Clustered Parameters**: R0, SI, CFR with uncertainty information
- **Full Distribution Parameters**: SI and IP distributions for presymptomatic calculation
- **Transmission Routes**: Binary indicators for each route type
- **Uncertainty Types**: Point estimates, confidence intervals, or ranges
- **Sampling Distributions**: Normal, beta, gamma, lognormal, or uniform

## Key Methodological Features

### Uncertainty Handling
- **Probabilistic Sampling**: Accounts for parameter uncertainty through MCMC
- **Multiple Distribution Types**: Supports various uncertainty representations
- **Robust Fallbacks**: Handles edge cases and missing data gracefully

### Clustering Approach
- **Ensemble Method**: Combines results across 5,000 MCMC iterations
- **Consensus Clustering**: Uses hierarchical clustering on co-assignment frequencies
- **Scaled Features**: Normalizes numerical parameters for clustering
- **Modal Assignments**: Determines most frequent cluster for each pathogen

### Validation Methods
- **Silhouette Analysis**: Evaluates clustering quality
- **Sensitivity Testing**: Assesses robustness to parameter inclusion/exclusion
- **Visual Assessment**: PCA plots and dendrograms for interpretation

## Output Interpretation

### Cluster Characteristics
- **Centroids**: Mean parameter values for each cluster
- **Confidence Intervals**: 95% CI across MCMC iterations
- **Route Profiles**: Average transmission route patterns

### Consensus Results
- **Stable Assignments**: Pathogens consistently grouped together
- **Cluster Stability**: Frequency of co-assignment across iterations
- **Hierarchical Structure**: Dendrogram showing pathogen relationships

## Usage Instructions

1. **Run Main Analysis**:
   ```r
   source("Clustering/mcmc/mcmc_Kmeans_analysis.R")
   ```

2. **Find Optimal K** (after main analysis):
   ```r
   source("Clustering/mcmc/find_optimal_K.R")
   ```

3. **Run Sensitivity Analysis**:
   ```r
   source("Clustering/mcmc/mcmc_sensitivity_analyis.R")
   ```

## Dependencies

Required R packages:
- `dplyr`, `readr`: Data manipulation
- `stats`: Statistical functions and clustering
- `cluster`: Silhouette analysis
- `ggplot2`, `factoextra`, `ggrepel`: Visualization
- `dendextend`, `RColorBrewer`: Enhanced dendrograms



## File Structure

```
Clustering/mcmc/
├── mcmc_Kmeans_analysis.R          # Main analysis script
├── find_optimal_K.R                # Optimal K determination
├── mcmc_sensitivity_analyis.R      # Sensitivity analysis
├── Kmeans/
│   ├── pathogen_params_*.csv       # Input data files
│   ├── mcmc_parameter_samples.csv  # MCMC results
│   ├── consensus_*.csv             # Clustering results
│   ├── *.png                       # Visualizations
│   └── MCMC_sens_analysis/         # Sensitivity analysis outputs
```

This project provides a robust, uncertainty-aware approach to pathogen classification that can inform public health decision-making and outbreak response strategies.
