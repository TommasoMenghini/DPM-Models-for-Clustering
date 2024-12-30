Description of the worldbank_data Application
================
As described in the [`README.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/README.md) file, this tutorial provides general guidelines and code to perform clustering using **Dirichlet Process Mixture (DPM) models** on the worldbank_data dataset.   
This file is specular to [`country_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/country_data.md) file, in fact it covers:

- How to **sample from a selection of multivariate Dirichlet Process Mixture models with Gaussian kernels**;
- Methods to **assess the convergence of the MCMC chain**;
- Implementation of the **SALSO algorithm with the VI loss function** to determine the optimal posterior partition of the data;
- Visualization and interpretation of **key graphical results** related to the clustering process;

Upload the worldbank Dataset
================

First set the working directory where `worldbank_data.csv` is placed. Once this has been done, clean the workspace, and load the data. The dataframe contains **socio-economic and health metrics** of 217 different nations. The first column lists the names of the different countries. The design matrix including the covariates can be easily obtained by extracting the remaining columns in `worldbank_data.csv`. In this case there are 7 covariates, two less than the `country_data.csv` example; it is not present the variable `Health`, that measures the percentage of per capita health expenditure, and `Income`, which represents the net per capita income. A possible direction for further study could be to understand how much the presence (or absence) of these two variables influences the categorization of a state and, consequently, the clustering process.

``` r


```
