Description of the Simulation
================
As described in the [`README.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/README.md) file, this tutorial provides general guidelines and code to perform clustering using **Dirichlet Process Mixture (DPM) models** on a simulated dataset. Specifically, it covers the following:

- Steps needed to create the **simulated data**;
- How to **sample from a selection of multivariate Dirichlet Process Mixture models with Gaussian kernels**;
- Methods to **assess the convergence of the MCMC chain**;
- Implementation of different **algorithms with VI and Binder loss functions** to determine the optimal posterior partition of the data;
- **Comparison** of different combinations of algorithms and loss functions to identify the one that provides the **most efficient solution for posterior inference**;

Generating the Simulated Dataset
================

In the simulation study we decided to generate 300 observation from a trivariate Gaussian mixture of 5 components. It can be written as:
$
a^2 + b^2 = c^2
$

``` r
rm(list=ls())

library(MASS)

set.seed(123)
n <- 300  
clusters <- 5 

mu_list <- list(
  c(-2, -2, -2), 
  c(2, -2, 2),
  c(-2, 2, -2),
  c(2, 2, 2),
  c(0,0,0)
)  

sigma_list <- list(
  diag(3),  
  diag(c(1.25, 1.25, 1.25)),
  diag(3),
  diag(c(1.25,1.25,1.25)),
  diag(c(0.5, 0.5, 0.5))
)

sim_data <- matrix(0, n, 4)

p <- rep(1/5, 5) 

for (i in 1:n) {
  # Selects one component of the mixture
  component <- sample(1:length(p), size = 1, prob = p)

  # Generates a sample from the distribution selected
  if (component == 1) {
    sim_data[i, 1:3] <- mvrnorm(1, mu = mu_list[[1]], Sigma = sigma_list[[1]])
    sim_data[i, 4] <- 1
  } else if (component == 2) {
    sim_data[i, 1:3] <- mvrnorm(1, mu = mu_list[[2]], Sigma = sigma_list[[2]])
    sim_data[i, 4] <- 2
  } else if (component == 3) {
    sim_data[i, 1:3] <- mvrnorm(1, mu = mu_list[[3]], Sigma = sigma_list[[3]])
    sim_data[i, 4] <- 3
  } else if(component == 4) {
    sim_data[i, 1:3] <- mvrnorm(1, mu = mu_list[[4]], Sigma = sigma_list[[4]])
    sim_data[i, 4] <- 4
  }  else{
    sim_data[i, 1:3] <- mvrnorm(1, mu = mu_list[[5]], Sigma = sigma_list[[5]])
    sim_data[i, 4] <- 5   
  }
}

sim_data <- as.data.frame(sim_data)
colnames(sim_data) <- c("X1", "X2", "X3", "Componente")
sim_data$Componente <- as.factor(sim_data$Componente)
x <- sim_data[,c(1:3)] 
cls.true <- sim_data[, 4]

```
