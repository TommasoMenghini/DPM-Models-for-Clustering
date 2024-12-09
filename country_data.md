Description of the Country_Data Application
================
As described in the [`README.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/README.md) file, this tutorial provides general guidelines and code to perform clustering using Dirichlet Process Mixture (DPM) models on the Country_Data dataset. Specifically, it covers the following:

- How to **sample from a selection of multivariate Dirichlet Process Mixture models with Gaussian kernels**;
- Methods to **assess the convergence of the MCMC chain**;
- Implementation of the **SALSO algorithm with the VI loss function** to determine the optimal posterior partition of the data;
- Visualization and interpretation of **key graphical results** related to the clustering process;

Upload the Country_Data Dataset
================

First set the working directory where `country-data.csv` is placed. Once this has been done, clean the workspace, and load the data. The dataframe contains socio-economic and health metrics of 167 different nations. The first column lists the names of the different countries. The design matrix including the covariates can be easily obtained by extracting the remaining columns in `country-data.csv`.

``` r
rm(list=ls())
data <- read.csv("country-data.csv")

country <- data[, 1]
x <- data[, -c(1)]

p <- dim(x)[2]
n <- dim(x)[1]

```

MCMC Sample Generation
================

Load the fundamental library `BNPmix` to implement the algorithm which provides i.i.d. samples from the latent partition posterior. First set the number of iterazions and burn-in of the MCMC chain. Then select the MCMC sampling method to be used, we were interested in the marginal sampler `MAR` by [`Neal 2000`](https://www.jstor.org/stable/1390653). The model is set to be `DLS` for a detailed and useful explanation see the [`BNPmix manual`](https://cran.r-project.org/web/packages/BNPmix). Lastly hyper is equal to FALSE so hyperprior distributions on the base measures's parameter are not added. 

Now set the hyperparameters of the that regulate the prior information. Setting the discount parameter to 0 is crucial, as it leads to a Dirichlet process, which is a specific case of the Pitman-Yor process for which the function `PYdensity()` is designed. The strength parameter is determined by a Gamma(1,1) random variable as suggested in [`Escobar and West 1995`](https://user-web-p-u02.wpi.edu/~balnan/Escobar-West-1995.pdf). The last argument is k0 that is the p-dimensional vector of scale factors defining the normal base measure on the location parameter and it is fixed empirically as a p-vector of 2. Matter of fact there are a lot of other parameters to modify, but we decided to leave the default values.

The last arguments to decide are those for generationg the posterior output. Being interested in the estimated partition put out_type equal to "CLUST".

``` r
library(BNPmix)

mcmc <- list(niter = 10000, nburn = 1000, method = "MAR", model = "DLS", 
             hyper = F)
prior <- list(strenght = rgamma(1,1,1), discount = 0, k0 = rep(2, p))
output <- list(out_type = 'CLUST', mean_dens = T)

set.seed(123)
fit <- PYdensity(y = x, mcmc = mcmc, prior = prior,
                 output = output)

```

MCMC Convergence Assessment
================

To assess convergence in our MCMC chain, we opted for a graphical approach, plotting two functionals of the chain: the number of clusters and the entropy of every visited partition. The resulting plots suggest convergence.

Another approach could involve applying diagnostics to these quantities. For instance, the R library `coda` offers a range of useful functions, such as `geweke.diag()`, which provides a convergence diagnostic based on a test for the equality of the means of the first and last parts of a Markov chain. If the test statistic is not significant, it is a positive indication of convergence.

``` r
library(coda)

partitions <- as.matrix(fit$clust)
num_clusters <- apply(partizioni, 1, function(partizioni) length(unique(partizioni)))

num_clusters <- as.mcmc(num_clusters)
plot(num_clusters)

geweke.diag(num_clusters)

H <- function(partition){
  N <- n
  kn <- max(partition)
  nj <- c()
  for (j in 1:kn){
    nj[j] <- sum(partition == j)
  }
  sum <- sum(nj*log(nj))
  log(N, base = 2) - 1/N*sum
}

entropy <- apply(partitions, 1, H)
entropy <- as.mcmc(entropy)
plot(entropy)

geweke.diag(entropy)

```

Salso + VI loss
================

As discussed in [`simulation.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/simulation.md), we compared different algorithms and loss functions to identify the combination that provides the most efficient solution for posterior inference. We found that the combination of the `SALSO` algorithm with the `VI` loss is the best choice for this problem. Therefore, we applied this approach to the real data.

``` r
library(salso)

start_time <- Sys.time()
salso.VI <- salso(fit$clust, loss=VI())
end_time <- Sys.time()

time_salso.VI <- difftime(end_time, start_time, units=("secs"))[[1]]

summ.VI <- summary(salso.VI)

labels.salso.VI <- summ.VI$estimate

pairs(x, col = labels.salso.VI, pch = 19)

```

Graphical Results
================

``` r

library(ggplot2)
library(rnaturalearthdata)
library(rnaturalearth)

labels.clust = as.numeric(summ.VI$estimate)
clust = cbind(cat, as.numeric(labels.clust))
colnames(clust) = c('country','cluster')
cluster = as.numeric(clust[,2])
clust = clust[,-2]
clust = as.data.frame(cbind(clust, cluster))
clust = clust[,c(1,2)]
clust[,2] = as.numeric(clust[,2])
colnames(clust) = c('country','cluster')

world = ne_countries(scale = "medium", returnclass = "sf")
unique(world$name)
unique(clust$country)

country_corrections = c(
  "Congo, Dem. Rep." = "Dem. Rep. Congo",
  "Congo, Rep." = "Congo",
  "Central African Republic" = "Central African Rep.",
  "United States" = "United States of America",
  "Czech Republic" =  "Czechia",
  "Slovak Republic" = "Slovakia",
  "Bosnia and Herzegovina" = "Bosnia and Herz.",
  "Bulgaria"="Bulgaria",
  "Kyrgyz Republic" = "Kyrgyzstan",
  "Cambodia" = "Cambodia",
  "Macedonia, FYR" = "North Macedonia",
  "St. Vincent and the Grenadines" = "St. Vin. and Gren.",
  "Equatorial Guinea"="Eq. Guinea",
  "Cote d'Ivoire" = "CÃ´te d'Ivoire"
)

clust = clust %>%
  mutate(country = recode(country, !!!country_corrections))

world = left_join(world, clust, by = c("name" = "country"))

ggplot(data = world) +
  geom_sf(aes(fill = factor(cluster))) +  
  scale_fill_manual(values = c("1" = 'yellow', "2" = 'red', "3" = 'green',
                               "4" = 'orange', "5" = 'purple', "6" = 'blue', na.value = 'grey'),
                    name = 'Cluster VI salso',
                    labels = c('1', '2', '3', '4', '5', '6', 'Valori non registrati')
  ) + 
  theme_minimal() +
  theme(legend.position = "bottom")


```

immagine grafico sopra


``` r

data.new <- as.data.frame(cbind(data, as.factor(clust$cluster)))
colnames(data.new)[colnames(data.new) == "as.factor(clust$cluster)"] <- "cluster"

scaled.gdp <- scale(data.new$gdpp)
df <- as.data.frame(cbind(as.numeric(scaled.gdp), as.factor(cluster), data$country))
colnames(df) = c('scaled.gdpp','cluster', "Country")

str(df)


which(data$country == "China")
which(data$country == "Chile")
which(data$country == "Ireland")
which(data$country == "Italy")
which(data$country == "Saudi Arabia")

      
density_estimate <- density(scaled.gdp)
# Interpola i valori di densità per il tuo vettore x
densities <- approx(density_estimate$x, density_estimate$y, xout = scaled.gdp)$y
length(densities)
length(data$country)

ggplot(df, aes(x = as.numeric(scaled.gdpp))) +
  geom_density() +  # Stima della densità
  geom_point(aes(y = densities , color = as.factor(cluster)),  # Punti colorati per cluster
             position = position_jitter(height = 0.01), size = 2) +
  scale_color_manual(values = c("1" = 'yellow', "2" = 'red', "3" = 'green',
                                "4" = 'orange', "5" = 'purple', "6" = 'blue', na.value = 'grey'),
                     name = 'Cluster VI salso',
                     labels = c('1', '2', '3', '4', '5', '6', 'Valori non registrati')
  ) +  # Usa il vettore di colori
  labs(x = "GDP pro capite scalato", y = "") +
  theme_minimal() +
  theme(legend.position = "none")  +
  # Aggiungi l'annotazione per la Cina
  annotate("text", x = 1+0.04, 
           y = densities[35] + 0.04, 
           label = "China", color = "black", size = 3, hjust = 0, fontface = "bold") +
  annotate("segment", x = as.numeric(df$scaled.gdpp[df$Country == "China"]), 
           xend = 1, 
           y = densities[35], 
           yend = densities[35]+0.02, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  annotate("text", x = 1+0.04, 
           y = 0.3+0.04, 
           label = "Saudi Arabia", color = "black", size = 3, hjust = 0, fontface = "bold") +
  annotate("segment", x = as.numeric(df$scaled.gdpp[df$Country == "Saudi Arabia"]), 
           xend = 1, 
           y = densities[129], 
           yend = 0.3, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  annotate("text", x = 3+0.04, 
           y = 0.3+0.04, 
           label = "Ireland", color = "black", size = 3, hjust = 0, fontface = "bold") +
  annotate("segment", x = as.numeric(df$scaled.gdpp[df$Country == "Ireland"]), 
           xend = 3, 
           y = densities[74], 
           yend = 0.3, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") 


labels.2010 <- clust$cluster
country2010 <- clust$country 

save(labels.2010, country2010, file = "2010data")

```
immagine grafico 





