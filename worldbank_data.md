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
Let's consider the graph below. The Dirichlet process mixture model identifies **8 different clusters**. When comparing these results with those reported in the [`country_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/country_data.md) file; it can be observed that the `Red` and `Yellow` clusters retain a similar interpretation. However, there are also significant shifts that highlight how the global landscape has changed over the twelve years separating this graph from the earlier one. **India** is a stricking example: in the 2010 graph it was classified as a `Yellow` country, but in the 2022 graph it had shifted to a `Red` one. his transition aligns with expectations, as India is a developing country in several respects—not only economically, but also culturally and politically. 

Before delving into further considerations similar to the previous one, it is necessary to make a brief disclaimer. There are no theoretical reasons to assume that the clustering obtained in the [`country_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/country_data.md) example should have any direct relationship with the clustering showed in the current graph. Even though this is undoubtedly true, it is also evident that some clusters share common countries across both examples, and therefore reflect common traits. For insance the `Red` cluster contains in both examples the countries that could be associated with the **Second Wold**, while the `Yellow` cluster is made up of poorer nations typically classified as part of the **Third World**. Encouraged by these parallels, we now turn our attention to three countries of particular interest: **China**, **Saudi Arabia** and **Ireland**. 

![](https://raw.githubusercontent.com/TommasoMenghini/DPM-Models-for-Clustering/main/img/World2022.png)

Run the following code to build a frequency table of the clusters.

``` r
library(knitr)

table(clust$cluster)[1]

tab <- c(table(clust$cluster)[1], table(clust$cluster)[2],
         table(clust$cluster)[3], table(clust$cluster)[4],
         table(clust$cluster)[5], table(clust$cluster)[6],
         table(clust$cluster)[6], table(clust

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

| Yellow| Red| Green| Orange| Purple| Blue| Pink| Violet|
|------:|---:|-----:|------:|------:|----:|----:|------:|
|     41|  47|    58|      5|      1|    4|    1|      2|


The resulting table shows the frequency of each cluster.

Scaled GDP per capita density plot
------------------

As introduced in the [`country_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/country_data.md) file let's now discuss about the **density of the scaled GDP per capita** plot, where each point represents a specific country and is colored according to the **2010 clustering partion**. The density of this feature was estimated using the R function `density()`, which computes kernel density estimates. Although alternative density estimation methods might be more elegant, we opted for this approach due to time constraints. As said before we address the case of  **China**, **Saudi Arabia** and **Ireland**. 

First we load some relevant quantities obtained from the [`country_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/country_data.md) file.

``` r

load("2010data.RData")

country2022 <- clust$country
country2010 == country2022
length(country2022); length(country2010)

```

Now we have to assign to each country in the 2022 dataset the label that was obtained in the 2010 study. The problem is that not all countries present in 2010 are also present in 2022 but first, we need to sort the countries in both the 2022 and 2010 datasets in alphabetical order.
This will be useful for correctly matching each label to the corresponding country.

``` r

country2022 <- sort(country2022)

dt <- data.frame(country2010, labels.2010)
dt <- dt %>% arrange(country2010)

index.common_nations <- which(dt$country2010 %in% country2022)
labels.2010_common <- dt$labels.2010[index.common_nations]
country2010_common <- dt$country2010[index.common_nations]
length(country2010_common)

```

Of course we need to consider only the GDP per capita of nations detected as common between the two years. Then it is straightforward to obtain the density estimation with the R function `density()`.

``` r

data$country.corr <- clust$country
dt <- data %>% select(gdpp, country.corr) %>% arrange(country.corr)
dt <- dt %>% filter(country.corr %in% country2010_common)

scaled.gdp <- scale(dt$gdpp)
length(scaled.gdp)

df <- as.data.frame(cbind(scaled.gdp, country2010_common, labels.2010_common))
colnames(df) <- c('scaled.gdpp','country', "cluster")

str(df)
df$scaled.gdpp <- as.numeric(df$scaled.gdpp)

which(df$country == "China")
which(df$country == "Ireland")
which(df$country == "Saudi Arabia")



density_estimate <- density(df$scaled.gdpp)
densities <- approx(density_estimate$x, density_estimate$y, xout = df$scaled.gdpp)$y

```

Finally we plot our results. The aim of this graph is to assess whether a country has shifted towards the right tail of the distribution, comparing its position in the plot below with the one on the corresponding plot in [`country_data.md`](https://github.com/TommasoMenghini/DPM-Models-for-Clustering/blob/main/country_data.md). Such a shift would suggest that the country has improved its economic, social and political conditions over the time interval considered. Indeed the wealthiest countries are situated closer to the right tail of the distributions. 

We can observe that **China** has made a significant shift toward the right tail of the distribution. This result is consistent with what is shown on the **World map** above: now **China** is a `Green` state, whereas in **2010** was labeled as a `Red` one. This outcome is not surprising as it is known that **China** undergone tremendous development in its socio-economic infrastructure and has assumed a central role in international politics.


``` r

ggplot(df, aes(x = scaled.gdpp)) +
  geom_density(size = 1.5) +  
  geom_point(aes(y = densities , color = as.factor(cluster)),  
             position = position_jitter(height = 0.01), size = 4) +
  scale_color_manual(values = c("1" = 'yellow', "2" = 'red', "3" = 'green',
                                "4" = 'orange', "5" = 'purple', "6" = 'blue', na.value = 'grey'),
                     name = 'Cluster VI salso',
                     labels = c('1', '2', '3', '4', '5', '6', 'No data available')
  ) +  
  labs(x = "Scaled GDP per capita",
       y = NULL) +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 20))  +
  annotate("text", x = df$scaled.gdpp[31]+0.55, 
           y = densities[31] + 0.04, 
           label = "China", color = "black", size = 10, hjust = 0, fontface = "bold") +
  annotate("segment", x = df$scaled.gdpp[31], 
           xend = df$scaled.gdpp[31]+0.5, 
           y = densities[31], 
           yend = densities[31]+0.02, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 2) +
  annotate("text", x = df$scaled.gdpp[110]+0.55, 
           y = densities[110]+0.04, 
           label = "Saudi Arabia", color = "black", size = 10, hjust = 0, fontface = "bold") +
  annotate("segment", x = df$scaled.gdpp[110], 
           xend = df$scaled.gdpp[110]+0.5, 
           y = densities[110], 
           yend = densities[110] + 0.02, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 2) +
  annotate("text", x = df$scaled.gdpp[65], 
           y =  densities[65]+0.25, 
           label = "Ireland", color = "black", size = 10, hjust = 0, fontface = "bold") +
  annotate("segment", x = df$scaled.gdpp[65], 
           xend = df$scaled.gdpp[65], 
           y = densities[65], 
           yend =  densities[65]+0.2, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 2) 

```



![](https://raw.githubusercontent.com/TommasoMenghini/DPM-Models-for-Clustering/main/img/GDPP_2022.png)










