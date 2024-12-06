Description of the Country_Data Application
================
As described in the README.md file, this tutorial provides general guidelines and code to perform clustering using Dirichlet Process Mixture (DPM) models on the Country_Data dataset. Specifically, it covers the following:

- How to **sample from a selection of multivariate Dirichlet Process Mixture models with Gaussian kernels**;
- Methods to **assess the convergence of the MCMC chain**;
- Implementation of the **SALSO algorithm with the VI loss function** to determine the optimal posterior partition of the data;
- Visualization and interpretation of **key graphical results** related to the clustering process;

Upload the Country_Data Dataset
================

breve descrizione

``` r
data <- read.csv("country-data.csv")
str(data)

x <- data[, -c(1)]
cat <- as.factor(data[, 1])

p <- dim(x)[2]
n <- dim(x)[1]

```

MCMC Sample Generation
================

Descrivo cosa faccio

``` r

library(BNPmix)
library(coda)

prior <- list(strenght = rgamma(1,1,1), discount = 0, k0 = rep(2, p))
output <- list(out_type = 'CLUST', mean_dens = T)
mcmc <- list(niter = 10000, nburn = 1000, method = "MAR", model = "DLS", 
             hyper = F)
set.seed(123)
fit <- PYdensity(y = x, mcmc = mcmc, prior = prior,
                 output = output)

```

MCMC Convergence Assessment
================

descrivo 

``` r
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

cosa mi risulta

Salso + VI loss
================

Descrivo cosa andrò a fare

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





