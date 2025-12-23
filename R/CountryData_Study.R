rm(list = ls())

data <- read.csv("country-data.csv")
str(data)


x <- data[, -c(1)]
country <- data[, 1]

p <- dim(x)[2]
n <- dim(x)[1]

# MCMC SAMPLE -----------------------------------------------

library(BNPmix)
library(coda)

prior <- list(strenght = rgamma(1,1,1), discount = 0, k0 = rep(2, p))
output <- list(out_type = 'CLUST', mean_dens = T)
mcmc <- list(niter = 10000, nburn = 1000, method = "MAR", model = "DLS", 
             hyper = F)
set.seed(123)
fit <- PYdensity(y = x, mcmc = mcmc, prior = prior,
                 output = output)

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



# SALSO + loss VI ---------------------------------------------------------

library(salso)

start_time <- Sys.time()
salso.VI <- salso(fit$clust, loss=VI())
end_time <- Sys.time()

time_salso.VI <- difftime(end_time, start_time, units=("secs"))[[1]]
# 1.866496

summ.VI <- summary(salso.VI)

labels.salso.VI <- summ.VI$estimate

pairs(x, col = labels.salso.VI, pch = 19)



# Graphical Results -------------------------------------------------------

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

country_corrections <- c(
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
  "Cote d'Ivoire" = "Côte d'Ivoire"
)

clust <- clust %>%
  mutate(country = recode(country, !!!country_corrections))

world <- left_join(world, clust, by = c("name" = "country"))

ggplot(data = world) +
  geom_sf(aes(fill = factor(cluster))) +  
  scale_fill_manual(values = c("1" = 'yellow', "2" = 'red', "3" = 'green',
                               "4" = 'orange', "5" = 'purple', "6" = 'blue', na.value = 'grey'),
                    name = 'Cluster VI salso',
                    labels = c('1', '2', '3', '4', '5', '6', 'No data available')
  ) + 
  theme_minimal() +
  theme(legend.position = "bottom")


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



table(clust$cluster)[1]
#  1  2  3  4  5  6 
# 58 66 29  7  4  3 
cluster6.2010 <- clust[clust$cluster == 6, 1]
# "Luxembourg" "Malta"      "Singapore" 
cluster5.2010 <- clust[clust$cluster == 5, 1]
# "Eq. Guinea" "Mongolia"   "Nigeria"    "Venezuela" 
cluster4.2010 <- clust[clust$cluster == 4, 1]



scaled.gdp <- scale(data$gdpp)
df <- as.data.frame(cbind(scaled.gdp, clust$cluster, data$country))
colnames(df) <- c('scaled.gdpp','cluster', "Country")

str(df)


which(data$country == "China")
which(data$country == "Chile")
which(data$country == "Ireland")
which(data$country == "Italy")
which(data$country == "Saudi Arabia")


density_estimate <- density(scaled.gdp)

densities <- approx(density_estimate$x, density_estimate$y, xout = scaled.gdp)$y
length(densities)
length(data$country)

ggplot(df, aes(x = as.numeric(scaled.gdpp))) +
  geom_density(size = 1.5) +  
  geom_point(aes(y = densities , color = as.factor(cluster)),  
             position = position_jitter(height = 0.01), size = 4) +
  scale_color_manual(values = c("1" = 'yellow', "2" = 'red', "3" = 'green',
                                "4" = 'orange', "5" = 'purple', "6" = 'blue', na.value = 'grey'),
                     name = 'Cluster VI salso',
                     labels = c('1', '2', '3', '4', '5', '6', 'No data available')
  ) +  
  labs(x = "Scaled GDP per capita", y = "") +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 20) )  +
  annotate("text", x = 1+0.04, 
           y = densities[35] + 0.04, 
           label = "China", color = "black", size = 10, hjust = 0, fontface = "bold") +
  annotate("segment", x = as.numeric(df$scaled.gdpp[df$Country == "China"]), 
           xend = 1, 
           y = densities[35], 
           yend = densities[35]+0.02, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 2) +
  annotate("text", x = 1+0.04, 
           y = 0.3+0.04, 
           label = "Saudi Arabia", color = "black", size = 10, hjust = 0, fontface = "bold") +
  annotate("segment", x = as.numeric(df$scaled.gdpp[df$Country == "Saudi Arabia"]), 
           xend = 1, 
           y = densities[129], 
           yend = 0.3, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 2) +
  annotate("text", x = 3+0.04, 
           y = 0.3+0.04, 
           label = "Ireland", color = "black", size = 10, hjust = 0, fontface = "bold") +
  annotate("segment", x = as.numeric(df$scaled.gdpp[df$Country == "Ireland"]), 
           xend = 3, 
           y = densities[74], 
           yend = 0.3, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 2) 

scaled_gdp.2010 <- scaled.gdp
labels.2010 <- clust$cluster
country2010 <- clust$country 

save(scaled_gdp.2010, labels.2010, country2010, file = "2010data.RData")

