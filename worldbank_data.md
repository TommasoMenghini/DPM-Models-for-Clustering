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

The last arguments to choose are those for generating the posterior output. Being **interested in the estimated partition** put out_type equal to "CLUST".

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

Graphical Results
================

This section contains the code to build two plots: 
- the **World Map**, that clearly shows the partition of countries that we obtained minimizing the posterior expected loss;
- the **Scaled GDP per Capita Density Plot**, which shows the estimated density of a scaled variable from the dataset. The optimal clustering estimate is also highlighted, with colours representing cluster membership.
  
World Map
------------------

Obtained the partition that minimizes the posterior expected loss, load some useful `R` packages in order to build the **World map**. Of course the names of certain countries need to be arranged. The following code generates the plot shown below.   

``` r
library(ggplot2)
library(rnaturalearthdata)
library(rnaturalearth)
library(tidyverse)

cat <- data[, 1]
labels.clust <- summ.VI$estimate

clust <- as.data.frame(cbind(cat, labels.clust))
clust[,2] <- as.numeric(clust[,2])
colnames(clust) <- c('country','cluster')
str(clust)

world <- ne_countries(scale = "medium", returnclass = "sf")
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
  "Cote d'Ivoire" = "CÃ´te d'Ivoire",
  "Russian Federation" = "Russia",
  "Egypt, Arab Rep." = "Egypt",
  "Niger" = "Nigeria",
  "Iran, Islamic Rep." = "Iran",
  "Turkiye" = "Turkey",
  "Viet Nam" = "Vietnam",
  "Korea, Rep." = "South Korea"
)

clust <- clust %>%
  mutate(country = recode(country, !!!country_corrections))

world <- left_join(world, clust, by = c("name" = "country"))

ggplot(data = world) +
  geom_sf(aes(fill = factor(cluster))) +  
  scale_fill_manual(values = c("1" = 'yellow', "2" = 'green', "3" = 'red',
                               "4" = 'orange', "5" = 'purple', "6" = 'blue',
                               "7" = 'pink', "8" = 'violet', na.value = 'grey'),
                    name = 'Cluster VI salso',
                    labels = c('1', '2', '3', '4', '5', '6', '7','8','Valori non registrati')
  ) + 
  theme_minimal() +
  labs(fill = "Cluster", title = "Anno 2022") +
  theme(legend.position = "bottom")
```
Let's consider the graph below. The Dirichlet process mixture model identifies **6 different clusters** and it's quite easy to recognize some patterns in the partition. Focus on the `Green` cluster, which includes nations that belong to the **West**, not in a geographical sense, but rather in a political and economic context. For example **Japan, South Corea and Australia** are part of this cluster but they are not located in Europe or North America, therefore we can say that these `Green` countries are those that are more developed in socio-economic and health terms.

The `Red` cluster contains countries that, at least in 2010, could be considered part of the **Second World** . These nations exhibit **varying levels of economic development**, often characterized by **strong industrialization** but also by **shortages in consumer goods** and living standards that are not consistently high. For example, the **Russian Federation** and **Brazil** are part of this cluster, which is convincing because these nations are known for economies heavily reliant on exporting raw materials such as natural gas and wood.

Focusing on the `Yellow` cluster, notice that it is spread across the Sub-Saharan region and the Indian subcontinent. Therefore, these states can be labeled as **Third World countries** —less developed economically, politically, and socially.

Let's talk about the `Orange` cluster: it contains only seven states — **Bahrain, Brunei, Kuwait, Oman, Qatar, Saudi Arabia, and the United Arab Emirates**. All these countries belong to the **Arabian Peninsula**, and they lead the global export of **oil**. That's why they are extremely wealthy nations. This wealth, however, is distributed in a highly unequal way: very few people, mainly the political and social elite, possess it. This is the most convincing argument that distinguishes these states from the `Green` ones. Furthermore, note that the only country in the Arabian Peninsula that does not belong to the `Orange` cluster is Yemen. Although **Yemen** possesses significant oil reserves, in 2010 it was not (and still is not) a major exporter of petroleum. The causes can be traced to internal difficulties such as civil war, corruption, and lacking infrastructure.

The `Purple` cluster contains four countries: Equatorial Guinea, Mongolia, Nigeria, and Venezuela. This time, there is no recognizable pattern.
The last cluster is the `Blue` one, which is not visible on the world map and contains three states: **Luxembourg, Malta, and Singapore**. In this case, it is fairly simple to understand what these nations have in common: they are geographically **very small** and can be (and still are) labeled as **tax havens**.


![](https://raw.githubusercontent.com/TommasoMenghini/DPM-Models-for-Clustering/main/img/World2010.png)

Run the following code to build a frequency table of the clusters.

``` r
library(knitr)

table(clust$cluster)[1]

tab <- c(table(clust$cluster)[1], table(clust$cluster)[2],
         table(clust$cluster)[3], table(clust$cluster)[4],
         table(clust$cluster)[5], table(clust$cluster)[6])

table <- data.frame(
  Yellow = tab[1],
  Red = tab[2],
  Green = tab[3],
  Orange = tab[4],
  Purple = tab[5],
  Blue = tab[6]
)

kable(table, format = "markdown")
```

| Yellow| Red| Green| Orange| Purple| Blue|
|------:|---:|-----:|------:|------:|----:|
|     58|  66|    29|      7|      4|    3|


The resulting table shows the frequency of each cluster, note that the largest clusters are the `Yellow` and `Red` ones.

Scaled GDP per capita density plot
------------------

In the [`worldbank_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/worldbakn_data.md) file, another World Map is created using data from 2022 instead of 2010, as shown in this example. Naturally, a different optimal partition that minimizes the posterior expected loss will be obtained. However, while it is not something to take for granted and caution is advised, **it is reasonable to compare clusters from different optimal partitions—the 2010 and the 2022 partitions**.

To pursue this objective, the plot below was created. It depicts the **density of the scaled GDP per capita**, where each point represents a specific country and is colored according to cluster membership. The density of this feature was estimated using the R function `density()`, which computes kernel density estimates. Although alternative density estimation methods might be more elegant, we opted for this approach due to time constraints.

Convinced by these arguments, we will study the cases of **China, Saudi Arabia, and Ireland**.

``` r

scaled.gdp <- scale(data$gdpp)
df <- as.data.frame(cbind(scaled.gdp, clust$cluster, data$country))
colnames(df) <- c('scaled.gdpp','cluster', "Country")

str(df)

which(data$country == "China")
which(data$country == "Ireland")
which(data$country == "Saudi Arabia")

      
density_estimate <- density(scaled.gdp)

densities <- approx(density_estimate$x, density_estimate$y, xout = scaled.gdp)$y
length(densities)
length(data$country)

ggplot(df, aes(x = as.numeric(scaled.gdpp))) +
  geom_density() +  
  geom_point(aes(y = densities , color = as.factor(cluster)),  
             position = position_jitter(height = 0.01), size = 2) +
  scale_color_manual(values = c("1" = 'yellow', "2" = 'red', "3" = 'green',
                                "4" = 'orange', "5" = 'purple', "6" = 'blue', na.value = 'grey'),
                     name = 'Cluster VI salso',
                     labels = c('1', '2', '3', '4', '5', '6', 'No data available')
  ) +  
  labs(x = "Scaled GDP per capita", y = "") +
  theme_minimal() +
  theme(legend.position = "none")  +
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

```

The inferences drawn from this plot will be discussed in the [`worldbank_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/worldbank_data.md) file. However, the main idea is to **observe the presence (or absence) of a shift to the right tail** of the distribution for a certain country, which could indicate its socio-economic development.


![](https://raw.githubusercontent.com/TommasoMenghini/DPM-Models-for-Clustering/main/img/GDPP_2010.png)










