Description of the Simulation
================
As described in the [`README.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/README.md) file, this tutorial provides general guidelines and code to perform clustering using **Dirichlet Process Mixture (DPM) models** on a simulated dataset. Specifically, it covers the following:

- Steps needed to create the simulated data;
- How to **sample from a selection of multivariate Dirichlet Process Mixture models with Gaussian kernels**;
- Methods to **assess the convergence of the MCMC chain**;
- Implementation of different **algorithms with VI and Binder loss functions** to determine the optimal posterior partition of the data;
- Visualization and interpretation of **key graphical results** related to the clustering process;

Upload the Country_Data Dataset
================

First set the working directory where `country-data.csv` is placed. Once this has been done, clean the workspace, and load the data. The dataframe contains **socio-economic and health metrics** of 167 different nations. The first column lists the names of the different countries. The design matrix including the covariates can be easily obtained by extracting the remaining columns in `country-data.csv`.

``` r
rm(list=ls())
data <- read.csv("country-data.csv")

country <- data[, 1]
x <- data[, -c(1)]

p <- dim(x)[2]
n <- dim(x)[1]

```
