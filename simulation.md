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

In the simulation study we decided to generate 300 observations from a **trivariate Gaussian mixture** of **5 components**. In the following table the means and variances of the components are shown:

<div align="center">


|               | j  =   1 |  j  =  2 | j = 3    | j = 4    | j = 5    |
|---------------|----------|----------|----------|----------|----------|
| $\mu_1$       |    -2    |     2    |   -2     |     2    |      0   |
| $\mu_2$       |     -2   |     -2   |    2     |    2     |     0    |
| $\mu_3$       |      -2  |      2   |    -2    |    2     |     0    |
| $\sigma^2$    |     1.00 |  1.25    |   1.00   |    1.25  |      0.50|

</div>

Below we provide the code to obtain such a simulated data.

``` r
rm(list=ls())

library(MASS)

set.seed(123)
n <- 300  
clusters <- 5 

mu_list <- list(c(-2, -2, -2), c(2, -2, 2), c(-2, 2, -2), c(2, 2, 2), c(0,0,0))  

sigma_list <- list(diag(3), diag(c(1.25, 1.25, 1.25)), diag(3), diag(c(1.25,1.25,1.25)), diag(c(0.5, 0.5, 0.5)))

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

It is useful to visualize the simulated data. Since the sample space is tridimensional, we suggest using interactive 3d plots to better understand the data. We used some useful libraries to help us in that sense and the code is shown below.

``` r
library(scatterplot3d)
library(plotly)
library(rgl)

# Interactive 3D plot
plot3d(x[,1], x[,2], x[,3], col = sim_data$Componente, type = "s", size = 1, xlab = "X1", ylab = "X2", zlab = "X3")

# Static 3D Plot
scatterplot3d(sim_data$X1, sim_data$X2, sim_data$X3, color = as.numeric(sim_data$Componente), pch = 19, main = "Grafico 3D statico dei dati simulati", xlab = "X1", ylab = "X2", zlab = "X3")

# Alternative Interactive 3D Plot
fig <- plot_ly(sim_data, x = ~X1, y = ~X2, z = ~X3, color = ~Componente, colors = "Set2", type = "scatter3d", mode = "markers", marker = list(size = 3))
fig <- fig %>% layout(scene = list(xaxis = list(title = 'X1'), yaxis = list(title = 'X2'), zaxis = list(title = 'X3')), title = "3D Interactive Scatterplot")
fig
```

The figure shown refers to the last plot presented in the code.
![](https://raw.githubusercontent.com/TommasoMenghini/DPM-Models-for-Clustering/main/img/Scatterplot3d.png)

MCMC Sample Generation
================

Load the fundamental library `BNPmix` to implement the algorithm which provides **i.i.d. samples** from the latent partition posterior. First set the **number of iterations** and **burn-in** of the MCMC chain. Then select the MCMC sampling method to be used, we were interested in the marginal sampler MAR by [`Neal 2000`](https://www.jstor.org/stable/1390653). The model is set to be **DLS** for a detailed and useful explanation see the [`BNPmix manual`](https://cran.r-project.org/web/packages/BNPmix). Lastly hyper is equal to FALSE so hyperprior distributions on the base measures's parameter are not added.

Now, we set the **hyperparameters** that govern the **prior information**. Setting the **discount parameter** to 0 is crucial, as this results in a **Dirichlet process**— a special case of the Pitman–Yor process for which the function `PYdensity()` is designed. The **strength parameter** is drawn from a `Gamma(1,1)` distribution, as suggested by [`Escobar and West 1995`](https://user-web-p-u02.wpi.edu/~balnan/Escobar-West-1995.pdf).

The parameters `m0` and `k0` represent the **mean vector** and the **vector of scale factors** that define the normal base measure on the location parameter. These were set as a three-dimensional vector of zeros and a three-dimensional vector of 1/2, respectively.

Similarly, `a0` and `b0` —which are **vectors of shape and scale parameters** defining the base measure on the scale parameters—were set as a three-dimensional vector of twos and a three-dimensional vector of ones, respectively.

Finally, execute the code below.

``` r
library(BNPmix)

prior <- list(strenght = rgamma(1,1,1), discount = 0, m0 = c(0,0,0), k0 = c(1/2,1/2,1/2),
              a0 = c(2,2,2), b0 = c(1,1,1))
output <- list(out_type = 'CLUST')
mcmc <- list(niter = 10000, nburn = 1000, method = "MAR", model = "DLS", 
             hyper = F)
set.seed(123)
fit <- PYdensity(y = x, mcmc = mcmc, prior = prior,
                 output = output)
```
The convergence of the MCMC sampler was assessed by tracking the entropy of the sampled partitions. Entropy provides a compact summary of partition complexity, jointly accounting for both the number of clusters and their relative sizes.

The corresponding traceplot shows stable oscillations after burn-in, with no evident trends, suggesting that the chain has reached a stationary regime.

Convergence was further checked using the Geweke diagnostic applied to the entropy chain. The test did not highlight significant differences between the early and late portions of the chain, supporting satisfactory convergence of the sampler, consistently with what already assessed.

``` r
H <- function(partizione){
  N <- n
  kn <- max(partizione)
  nj <- c()
  for (j in 1:kn){
    nj[j] <- sum(partizione == j)
  }
  somma <- sum(nj*log(nj))
  log(N, base = 2) - 1/N*somma
}

entropy <- apply(partizioni, 1, H)
entropy <- as.mcmc(entropy)
plot(entropy)

geweke.diag(entropy)
```

Implementation of different Algorithms with Loss Functions
================

Once convergence of the MCMC sampler has been assessed, the posterior distribution of partitions is summarized by minimizing different loss functions.
In this study, we consider two loss functions — **Variation of Information (VI)** and **Binder loss** — combined with two optimization strategies: **Greedy Search** [`Wade and Ghahramani`](https://projecteuclid.org/journals/bayesian-analysis/volume-13/issue-2/Bayesian-Cluster-Analysis--Point-Estimation-and-Credible-Balls-with/10.1214/17-BA1073.pdf) and **SALSO** [`Dahl et al.`](https://www.tandfonline.com/doi/abs/10.1080/10618600.2022.2069779).

All methods are applied starting from the posterior similarity matrix (PSM), computed from the MCMC output.

``` r
library(mcclust.ext)
psm <- comp.psm(fit$clust)
```

## Greedy Search with Variation of Information

The Greedy algorithm searches locally on the lattice of partitions to minimize the posterior expected VI loss. The resulting partition identifies five clusters, which matches the true number of components used in the data–generating process, we can see them in the table right below.

``` r
wade.VI <- minVI(psm, fit$clust, method = "all", include.greedy = TRUE)
labels.wade.VI <- wade.VI$cl[5, ]
```

|                                | j  =   1 |  j  =  2 | j = 3    | j = 4    | j = 5    |
|--------------------------------|----------|----------|----------|----------|----------|
| Real Cluster                   |    63    |     65   |   54     |     74   |     44   |
| VI Greedy Search Cluster       |     56   |     66   |    58    |    72    |     48   |

## Greedy Search with Binder loss

The same Greedy strategy is applied using the Binder loss with equal weights for the two types of misclassification. In this case, the estimated partition contains a substantially larger number of clusters, highlighting the tendency of the Binder loss to over–partition the data when uncertainty is present near cluster boundaries.

``` r
wade.B <- minbinder.ext(psm, fit$clust, method = "all", include.greedy = TRUE)
labels.wade.B <- wade.B$cl[5, ]
```
## SALSO with Variation of Information

The SALSO algorithm provides a stochastic alternative to Greedy Search and is designed to reduce sensitivity to local minima.
It directly optimizes the posterior expected loss using multiple randomized runs. Under the VI loss, SALSO again recovers a partition with five clusters, showing strong agreement with the Greedy–VI solution but at a much lower computational cost.

``` r
library(salso)

salso.VI <- salso(fit$clust, loss = VI())
labels.salso.VI <- summary(salso.VI)$estimate
```
