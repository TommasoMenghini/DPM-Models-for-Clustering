# DPM-Models-for-Clustering

# DPM Models for Bayesian Clustering

This repository contains materials developed for the **Nonparametric Statistics** exam project, focused on **Bayesian clustering via Dirichlet Process Mixture (DPM) models**.  
The project investigates **decision–theoretic approaches** to summarize posterior distributions over partitions, with particular attention to **loss functions** and **search algorithms**.

The work combines **theoretical insights**, **simulation studies**, and **real–data applications**, with fully reproducible **R code**.

---

## Background and References

Our analysis builds primarily on the following works:

- [`Wade and Gharamani `](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/README.md)
- [`Dahl et al.`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/README.md)

The exploration focuses on two loss functions:

- **Binder loss** (and its generalized version)
- **Variation of Information (VI)** loss

and on two optimization strategies:

- **Greedy Search**
- **SALSO (Search Algorithm for Loss-based clustering Optimization)**

---

## Repository Structure

The repository is organized as follows:



### Core scripts

- **`SimulationStudyTRUE.R`**  
  R code for the simulation study based on synthetic data generated from a multivariate Gaussian mixture.

- **`WorldBankStudy.R`**  
  Application of DPM models to World Bank socio–economic data (year 2022).

- **`CountryData_Study.R`**  
  Application of DPM models to the Country-Data dataset (HELP International, year 2010).

### Documentation

- **`simulation.md`**  
  Detailed description of the simulation design, model specification, convergence diagnostics, loss functions, algorithms, and performance comparison.

- **`worldbank_study.md`**  
  Explanation of the clustering analysis on World Bank data, including model assumptions, posterior inference, and interpretation of clusters.

- **`countrydata_study.md`**  
  Description of the analysis on the Country-Data dataset and comparison with more recent data.

Each `.md` file alternates **theoretical explanations**, **methodological choices**, and **selected R code snippets**, following a tutorial-style layout.

---

## Methodological Overview

### Model

All analyses rely on **Dirichlet Process Mixture Models** with Gaussian kernels.  
Inference on latent partitions is performed via **Marginal Gibbs Sampling** (Neal, 2000), as implemented in the `BNPmix` R package.

### Posterior Summarization

Given the posterior distribution over partitions, we adopt a **decision–theoretic framework**:

where \( L \) is either:

- Binder loss  
- Variation of Information loss  

### Optimization Algorithms

The expected posterior loss is minimized using:

- **Greedy Search** (Wade & Ghahramani, 2018)
- **SALSO** (Dahl et al., 2022)

We also study **credible balls** to quantify uncertainty around the estimated optimal partition.

---

## Simulation Study

A simulation study is conducted using trivariate Gaussian mixtures with known cluster structure.  
The study compares:

- Binder vs VI loss  
- Greedy vs SALSO algorithms  
- Performance across increasing sample sizes  

Metrics include:

- Number of recovered clusters
- Misclassified observations
- VI and Binder loss w.r.t. the true partition
- Adjusted Rand Index (ARI)
- Computational time

Results show that **SALSO combined with VI loss** provides the most stable and scalable performance.

---

## Real Data Applications

### Country-Data (2010)

Clustering of 167 countries based on socio–economic and health indicators, with the goal of identifying groups of nations in greatest need of humanitarian aid.

### World Bank Data (2022)

Extension of the analysis to more recent data, allowing for a qualitative comparison of global development patterns over time.

Both analyses reveal **interpretable and coherent cluster structures**, consistent with economic and geopolitical considerations.

---

## Computational Details

All analyses were performed using:

- **R version 4.3.3**
- **64-bit Windows environment**

Main R packages used include:

- `BNPmix`
- `salso`
- `mcclust.ext`
- `mclust`
- `coda`
- `ggplot2`
- `rnaturalearth`

---

## Reproducibility

All results presented in the `.md` files can be fully reproduced by running the corresponding `.R` scripts.  
Random seeds are fixed where appropriate.

---

## Authors

- **Tommaso Menghini**
- **Tommaso Pozzi**

Project developed for the **Nonparametric Statistics** exam.

