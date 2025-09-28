# R_code_spides_ind_niche
In the folder, there are two R codes used in "Contrasting spatial signatures of thermal and trophic individual-level niches in a community of tropical spiders" ("T_test.R" and "Spiders_var.R") and an RDS.file

# README: Thermal and Trophic Niches Analysis

## Overview
This project analyzes the thermal and trophic niches of individual species, assessing spatial autocorrelation and applying Generalized Least Squares (GLS) models to account for spatial dependencies.

------

# File "T-test.R"

This file contains the resampling protocol used. The code is written in a way that creates an RDS file to allow reproduction of the previously obtained results. If the code is run with the RDS file in the directory, it will replicate the same sampling indices, ensuring that all resampled datasets, statistical tests, and conclusions remain identical across runs. If the RDS file is not present in the directory, it will generate new random sampling indices based on the specified parameters, save them as a new RDS file.

------

# File "Spiders_var.R"

This file contains the code used to obtained the results. More below.

## Testing the Variables
Because variance analysis assumes data independence, we assessed whether the trophic and thermal niches of individual species were spatially autocorrelated. We used Moran’s I statistic via the `ape` package to determine if individuals with similar niches are spatially clustered.

### Spatial Autocorrelation
We tested spatial autocorrelation for:
- **Average temperature** (mean of all recorded temperatures for each individual)
- **Temperature width** (difference between minimum and maximum temperature)
- **δ13C** and **δ15N** (stable isotope ratios)

A significant spatial autocorrelation was detected for average temperature and δ13C, leading us to incorporate proximity matrices into our GLS models.

## Thermic Variables
### Moran's I Test for Spatial Autocorrelation
1. Organize longitude (`long`) and latitude (`lat`) from the dataset.
2. Compute a distance matrix using Euclidean distance.
3. Create an inverse distance matrix (proximity matrix) with diagonal values set to zero.
4. Apply Moran’s I test on log-transformed isotope values.

**Results:**
- **Average Temperature:** Significant spatial autocorrelation detected.
- **Temperature Width:** No significant spatial autocorrelation detected.

### Generalized Least Squares (GLS) Modeling
To account for spatial dependencies, we fitted GLS models using different correlation structures (Exponential, Gaussian, Ratio, Spherical) and selected the best model based on AIC values.

**Results:**
- The Gaussian correlation model was the best fit for **Average Temperature**.
- No spatial structure was needed for **Temperature Width**.
- ANOVA results showed species identity, body size, and their interaction significantly affected thermal niches.
- Tukey-adjusted pairwise contrasts identified significant differences between species.

## Trophic Variables
### δ13C and δ15N Moran’s I Tests
We performed Moran’s I tests on the isotopic data:
- **δ13C:** Significant spatial autocorrelation detected.
- **δ15N:** No significant spatial autocorrelation detected.

### GLS Modeling
To address spatial structure in δ13C, we fitted GLS models, selecting the optimal correlation structure based on AIC values.

**Results:**
- Gaussian correlation model best fit for **δ13C**.
- No spatial structure required for **δ15N**.
- Species identity significantly influenced both isotopic values.
- Pairwise contrasts revealed species-specific differences in trophic niches.

## Visualizations
We created density ridge plots to visualize the temperature distribution over time for different species using `ggplot2`, `ggridges`, and `viridis`. Data were reshaped to long format before plotting.

## Mantel test of spatial and niche distance matrices

To explore the potential covariation between body temperature, stable isotope values, and spatial position, we conducted three different Mantel tests
(accounting for each pairwise combination of variables) using the R package ecodist (Goslee & Urban 2007). Three Euclidean distance matrices were
generated: a thermal matrix (based on individual temperature records), a trophic matrix (using δ13C and δ15N values for all individuals), and a
spatial matrix (derived from the geographic coordinates of all individuals). The statistical significance of all Mantel tests was obtained by comparing
the Euclidean distance matrices of each factor with 10,000 random permutations of the arrangement of the cells within these matrices, assessing significant
differences from zero through Spearman correlations.

