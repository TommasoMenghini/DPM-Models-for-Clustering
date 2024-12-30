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

First set the working directory where `worldbank_data.csv` is placed. Once this has been done, clean the workspace, and load the data. The dataframe contains **socio-economic and health metrics** of 159 different nations. The first column lists the names of the different countries. The design matrix including the covariates can be easily obtained by extracting the remaining columns in `worldbank_data.csv`. In this case there are 7 covariates, two less than the `country_data.csv` example; it is not present the variable `Health`, that measures the percentage of per capita health expenditure, and `Income`, which represents the net per capita income. A possible direction for further study could be to understand how much the presence (or absence) of these two variables influences the categorization of a state and, consequently, the clustering process.

It's important to handle the NAs. We decided to omitt all those countries with at least one NA, but you can think to a smarter way to manage them.

``` r

rm(list = ls())
data <- read.csv('worldbank_data.csv')


colnames(data) <- c('country', 'contry_name', 'gdp', 'child_mort', 'exports', 'import', 'inflation',
                   'life_expec', 'total_fer', 'gdpp')


data$child_mort <- as.numeric(data$child_mort)
data$exports <- as.numeric(data$exports)
data$import <- as.numeric(data$import)
data$inflation <-  as.numeric(data$inflation)
data$life_expec <- as.numeric(data$life_expec)
data$total_fer <- as.numeric(data$total_fer)
data$gdpp <- as.numeric(data$gdpp)


data_clean <- na.omit(data)
data_clean <- data_clean[1:159, ] 

data <- data_clean[,-c(2,3)]

x <- round(data[,-1], 1) 

p <- dim(x)[2]
n <- dim(x)[1]

```
