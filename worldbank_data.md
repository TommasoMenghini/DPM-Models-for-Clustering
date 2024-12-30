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

MCMC Sample Generation
================

Load the fundamental library `BNPmix` to implement the algorithm which provides **i.i.d. samples** from the latent partition posterior. First set the number of iterations and burn-in of the MCMC chain. Then select the MCMC sampling method to be used, we were interested in the **marginal sampler** `MAR` by [`Neal 2000`](https://www.jstor.org/stable/1390653). The model is set to be `DLS` for a detailed and useful explanation see the [`BNPmix manual`](https://cran.r-project.org/web/packages/BNPmix). Lastly hyper is equal to `FALSE` so hyperprior distributions on the base measures's parameter are not added. 

Now set the **hyperparameters** that regulate the prior information. Setting the **discount parameter** to 0 is crucial, as it leads to a Dirichlet process, which is a specific case of the Pitman-Yor processes for which the function `PYdensity()` is designed. The **strength parameter** is determined by a Gamma(1,1) random variable as suggested in [`Escobar and West 1995`](https://user-web-p-u02.wpi.edu/~balnan/Escobar-West-1995.pdf). The last argument is k0 that is the p-dimensional vector of **scale factors** defining the **normal base measure** on the **location parameter** and it is fixed empirically as a **p-vector of 1/2**. Actually there are a lot of other parameters to modify, but we decided to leave the default values.

The last arguments to decide are those for generating the posterior output. Being **interested in the estimated partition** put out_type equal to "CLUST".

``` r
library(BNPmix)

mcmc <- list(niter = 20000, nburn = 5000, method = "MAR", model = "DLS", 
             hyper = F)
prior <- list(strenght = rgamma(1,1,1), discount = 0, k0 = rep(1/2, p))
output <- list(out_type = 'CLUST')

set.seed(123)
fit <- PYdensity(y = x, mcmc = mcmc, prior = prior,
                 output = output)

```

MCMC Convergence Assessment
================

To **assess convergence** in our MCMC chain you can opt for a graphical approach, plotting two functionals of the chain: the **number of clusters** and the **entropy** of every visited partition. 

Another approach could involve applying **diagnostics** to these quantities. For instance, the R library `coda` offers a range of useful functions, such as `geweke.diag()`, which provides a **convergence diagnostic based on a test** for the equality of the means of the first and last parts of a Markov chain. If the test statistic is not significant, it is a positive indication of convergence.

``` r
library(coda)

partizioni <- as.matrix(fit$clust)
num_clusters <- apply(partizioni, 1, function(partizioni) length(unique(partizioni)))

num_clusters <- as.mcmc(num_clusters)
plot(num_clusters)

geweke.diag(num_clusters)


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
The trace plots and the tests suggest convergence.


Salso + VI loss
================

As discussed in [`simulation.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/simulation.md), we compared **different algorithms** and **loss functions** to identify the combination that provides the **most efficient solution for posterior inference**. We found out that the combination of the `SALSO` algorithm with the `VI` loss is the best choice for this problem. Therefore, we applied this approach to the dataset obtaining the partition that minimizes the posterior expexted loss.

``` r
library(salso)


start_time <- Sys.time()
salso.VI <- salso(fit$clust, loss = salso::VI())
end_time <- Sys.time()

time_salso.VI <- difftime(end_time, start_time, units=("secs"))[[1]]

summ.VI <- summary(salso.VI)

labels.salso.VI <- summ.VI$estimate

pairs(x, col = labels.salso.VI, pch = 19)

```











