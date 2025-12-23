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


# MCMC SAMPLE -----------------------------------------------

library(BNPmix)

prior <- list(strenght = rgamma(1,1,1), discount = 0, k0 = rep(1/2, p))
output <- list(out_type = 'CLUST')
mcmc <- list(niter = 20000, nburn = 5000, method = "MAR", model = "DLS", 
             hyper = F)
set.seed(123)
fit <- PYdensity(y = x, mcmc = mcmc, prior = prior,
                 output = output)

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


# SALSO + loss VI ---------------------------------------------------------


library(salso)


start_time <- Sys.time()
salso.VI <- salso(fit$clust, loss = salso::VI())
end_time <- Sys.time()

time_salso.VI <- difftime(end_time, start_time, units=("secs"))[[1]]

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
  "Cote d'Ivoire" = "Côte d'Ivoire",
  "Russian Federation" = "Russia",
  "Egypt, Arab Rep." = "Egypt",
  "Niger" = "Nigeria",
  "Iran, Islamic Rep." = "Iran",
  "Turkiye" = "Turkey",
  "Viet Nam" = "Vietnam",
  "Korea, Rep." = "South Korea"
)


clust = clust %>%
  mutate(country = recode(country, !!!country_corrections))

world = left_join(world, clust, by = c("name" = "country"))

ggplot(data = world) +
  geom_sf(aes(fill = factor(cluster))) +  
  scale_fill_manual(values = c("1" = 'yellow', "2" = 'green', "3" = 'red',
                               "4" = 'orange', "5" = 'purple', "6" = 'blue',
                               "7" = 'pink', "8" = 'violet', na.value = 'grey'),
                    name = 'Cluster VI salso',
                    labels = c('1', '2', '3', '4', '5', '6', '7','8','Not Available')
  ) + 
  theme_minimal() +
  labs(fill = "Cluster", title = "Year 2022") +
  theme(legend.position = "bottom")


library(knitr)

table(clust$cluster)[1]

tab <- c(table(clust$cluster)[1], table(clust$cluster)[2],
         table(clust$cluster)[3], table(clust$cluster)[4],
         table(clust$cluster)[5], table(clust$cluster)[6],
         table(clust$cluster)[7], table(clust$cluster)[8])
                                        
                                        
table <- data.frame(
  Yellow = tab[1],
  Red = tab[3],
  Green = tab[2],
  Orange = tab[4],
  Purple = tab[5],
  Blue = tab[6],
  Pink = tab[7],
  Violet = tab[8]
)

kable(table, format = "markdown")


load("2010data.RData")

country2022 <- clust$country
country2010 == country2022
length(country2022); length(country2010)

labels.2010

country2022 <- sort(country2022)

dt <- data.frame(country2010, labels.2010)
dt <- dt %>% arrange(country2010)

index.common_nations <- which(dt$country2010 %in% country2022)
labels.2010_common <- dt$labels.2010[index.common_nations]
country2010_common <- dt$country2010[index.common_nations]
length(country2010_common)

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

