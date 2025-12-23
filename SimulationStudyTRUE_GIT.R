rm(list = ls())

library(MASS)   
library(scatterplot3d)
library(plotly)
library(BNPmix)
library(salso)
library(mclust)
library(coda)
library(rgl)


# Simulated Data -----------------------------------------------

set.seed(123)
n <- 300  
clusters <- 5 

mu_list <- list(
  c(-2, -2, -2), 
  c(2, -2, 2),
  c(-2, 2, -2),
  c(2, 2, 2),
  c(0,0,0)) 

sigma_list <- list(
  diag(3),  
  diag(c(1.25, 1.25, 1.25)),
  diag(3),
  diag(c(1.25,1.25,1.25)),
  diag(c(0.5, 0.5, 0.5)))


sim_data <- matrix(0, n, 4)

p <- rep(1/5, 5) 

for (i in 1:n) {
  
  component <- sample(1:length(p), size = 1, prob = p)
  
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
table(cls.true)


plot3d(x[,1], x[,2], x[,3], col = sim_data$Componente,
       type = "s", size = 1, xlab = "X1", ylab = "X2", zlab = "X3")

scatterplot3d(sim_data$X1, sim_data$X2, sim_data$X3, color = as.numeric(sim_data$Componente), pch = 19,
              main = "Grafico 3D statico dei dati simulati", xlab = "X1", ylab = "X2", zlab = "X3")

fig <- plot_ly(sim_data, x = ~X1, y = ~X2, z = ~X3, color = ~Componente, colors = "Set2", 
               type = "scatter3d", mode = "markers", marker = list(size = 3))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'X1'),
                                   yaxis = list(title = 'X2'),
                                   zaxis = list(title = 'X3')),
                      title = "3D Interactive Scatterplot")
fig



# MCMC SAMPLE -----------------------------------------------

prior <- list(strenght = rgamma(1,1,1), discount = 0, m0 = c(0,0,0), k0 = c(1/2,1/2,1/2),
              a0 = c(2,2,2), b0 = c(1,1,1))
output <- list(out_type = 'CLUST')
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

save(fit, cls.true, n, file= "MCMCfit.RData")

# Greedy Wade + VI --------------------------------------------------------

load("WorkSpace/MCMCfit.RData")

library(mcclust.ext)
library(mclust)
library(salso)
library(ggplot2)
library(plotly)

psm <- comp.psm(fit$clust)

start_time <- Sys.time()
wade.VI <- minVI(psm,fit$clust, method=("all"), include.greedy = T)
end_time <- Sys.time()

time_wade.VI <- difftime(end_time, start_time, units=("secs"))[[1]]

summary(wade.VI)


labels.wade.minVI <- wade.VI$cl[5, ]
plot(wade.VI, data = x, pch = 19)

pairs(x, col = labels.wade.minVI, pch = 19)
pairs(x, col = cls.true, pch = 19)

table(labels.wade.minVI)
table(cls.true)

labels.wade.minVI_assigned <- c()

for (j in 1:length(labels.wade.minVI)){
  if (labels.wade.minVI[j] == 1) {
    labels.wade.minVI_assigned[j] <- 3
  } else if (labels.wade.minVI[j] == 2) {
    labels.wade.minVI_assigned[j] <- 1
  } else if (labels.wade.minVI[j] == 3) {
    labels.wade.minVI_assigned[j] <- 2
  } else if (labels.wade.minVI[j] == 4){
    labels.wade.minVI_assigned[j] <- 5
  } else {
    labels.wade.minVI_assigned[j] <- 4
  }
}

table(labels.wade.minVI_assigned)
table(sim_data$Componente)

labels.wade.minVI_assigned == cls.true
pairs(x, col = labels.wade.minVI_assigned, pch = 19)

wade.cb.VI <- credibleball(wade.VI$cl[5,],fit$clust, c.dist = "VI", alpha = 0.05)
summary(wade.cb.VI)
plot(wade.cb.VI, data = x, pch = 19)

loss.wade_VIVI <- partition.loss(cls.true, labels.wade.minVI, loss = VI())
loss.wade_VIB <- partition.loss(cls.true, labels.wade.minVI, loss = binder())

ari_VI <- adjustedRandIndex(cls.true, labels.wade.minVI)

x.assigned.VI <- data.frame(cbind(x, as.factor(labels.wade.minVI_assigned)))
colnames(x.assigned.VI) <- c("X1", "X2", "X3", "Componente.VI")


p12 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 1], y = x.assigned.VI[, 2], color = Componente.VI)) +
  geom_point(size = 2) +               
  scale_color_brewer(palette = "Set2") +  
  labs(x = "X1", y = "X2", color = "Cluster") +      
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank())  

p13 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 1], y = x.assigned.VI[, 3], color = Componente.VI)) +
  geom_point(size = 2) +               
  scale_color_brewer(palette = "Set2") +  
  labs(x = "X1", y = "X3", color = "Cluster") +      
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank())  

p21 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 2], y = x.assigned.VI[, 1], color = Componente.VI)) +
  geom_point(size = 2) +               
  scale_color_brewer(palette = "Set2") +  
  labs(x = "X2", y = "X1", color = "Cluster") +      
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank())  

p23 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 2], y = x.assigned.VI[, 3], color = Componente.VI)) +
  geom_point(size = 2) +              
  scale_color_brewer(palette = "Set2") +  
  labs(x = "X2", y = "X3", color = "Cluster") +      
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank())  

p31 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 3], y = x.assigned.VI[, 1], color = Componente.VI)) +
  geom_point(size = 2) +               
  scale_color_brewer(palette = "Set2") +  
  labs(x = "X3", y = "X1", color = "Cluster") +      
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank())  

p32 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 3], y = x.assigned.VI[, 2], color = Componente.VI)) +
  geom_point(size = 2) +              
  scale_color_brewer(palette = "Set2") +  
  labs(x = "X3", y = "X2", color = "Cluster") +      
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank())  

p1 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 1])) +
  geom_density(fill = "lightblue", color = "blue") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank(), axis.title.y = element_blank())  

p2 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 2])) +
  geom_density(fill = "lightblue", color = "blue") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_blank(), axis.title.y = element_blank())  

p3 <- ggplot(x.assigned.VI, aes(x = x.assigned.VI[, 3])) +
  geom_density(fill = "lightblue", color = "blue") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_blank(), axis.title.y = element_blank())  

gridExtra::grid.arrange(p1,p12,p13,p21,p2,p23,p31,p32,p3, layout_matrix = matrix(1:9, 3,3))


sim_data.VI <- cbind(sim_data[, -4], labels.wade.minVI)
fig <- plot_ly(sim_data, x = ~X1, y = ~X2, z = ~X3, color = ~labels.wade.minVI, colors = "Set2", type = "scatter3d", mode = "markers", marker = list(size = 3))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'X1'),
                                   yaxis = list(title = 'X2'),
                                   zaxis = list(title = 'X3')),
                      title = "Grafico 3D interattivo dei dati simulati")
fig


save(labels.wade.minVI_assigned, loss.wade_VIVI, loss.wade_VIB, ari_VI, time_wade.VI, file= "Greedy_WadeVI.RData")


# Greedy Wade + Loss Binder  ---------------------------------

load("MCMCfit.RData")

library(mcclust.ext)
library(salso)
library(mclust)

psm <- comp.psm(fit$clust)
start_time <- Sys.time()
wade.B  <- minbinder.ext(psm, fit$clust, method=("all"), include.greedy=TRUE)
end_time <- Sys.time()

time_wade.B <- difftime(end_time, start_time, units=("secs"))[[1]]

summary(wade.B)
plot(wade.B, data = x)


labels.wade.minbinder <- wade.B$cl[5, ]

table(labels.wade.minbinder)
table(cls.true)

pairs(x, col = labels.wade.minbinder, pch = 19)
pairs(x, col = cls.true, pch = 19)


labels.wade.minbinder_assigned <- labels.wade.minbinder

for (j in 1:length(labels.wade.minbinder)){
  if (labels.wade.minbinder[j] == 1) {
    labels.wade.minbinder_assigned[j] <- 3
  } else if (labels.wade.minbinder[j] == 2) {
    labels.wade.minbinder_assigned[j] <- 1
  } else if (labels.wade.minbinder[j] == 3) {
    labels.wade.minbinder_assigned[j] <- 2
  } else if (labels.wade.minbinder[j] == 4) {
    labels.wade.minbinder_assigned[j] <- 5
  } else if (labels.wade.minbinder[j] == 5){
    labels.wade.minbinder_assigned[j] <- 4
  } 
}

table(labels.wade.minbinder_assigned)
table(cls.true)

pairs(x, col = cls.true, pch = 19)
pairs(x, col = labels.wade.minbinder_assigned, pch = 19)


wade.cb.B <- credibleball(wade.B$cl[5,], fit$clust)
summary(wade.cb.B)
plot(wade.cb.B, data = x, pch = 19)

labels.wade.minbinder <- wade.B$cl[5, ]

loss.wade_BVI <- partition.loss(cls.true, labels.wade.minbinder, loss = VI())
loss.wade_BB <- partition.loss(cls.true, labels.wade.minbinder, loss = binder())

ari_B <- adjustedRandIndex(cls.true, labels.wade.minbinder)

save(labels.wade.minbinder_assigned, loss.wade_BVI, loss.wade_BB, ari_B, time_wade.B, file= "Greedy_WadeB.RData")


# SALSO + loss VI ---------------------------------------------------------

load("MCMCfit.RData")

library(salso)
library(ggplot2)
library(mclust)

start_time <- Sys.time()
salso.VI <- salso(fit$clust, loss=VI())
end_time <- Sys.time()

time_salso.VI <- difftime(end_time, start_time, units=("secs"))[[1]]


summ.VI <- summary(salso.VI)

labels.salso.VI <- summ.VI$estimate

table(labels.salso.VI)
table(cls.true)

pairs(x, col = labels.salso.VI, pch = 19)
pairs(x, col = cls.true, pch = 19)


labels.salso.VI_assigned <- labels.salso.VI

for (j in 1:length(labels.salso.VI)){
  if (labels.salso.VI[j] == 1) {
    labels.salso.VI_assigned[j] <- 3
  } else if (labels.salso.VI[j] == 2) {
    labels.salso.VI_assigned[j] <- 1
  } else if (labels.salso.VI[j] == 3) {
    labels.salso.VI_assigned[j] <- 2
  } else if (labels.salso.VI[j] == 4) {
    labels.salso.VI_assigned[j] <- 5
  } else if (labels.salso.VI[j] == 5){
    labels.salso.VI_assigned[j] <- 4
  } 
}

table(labels.salso.VI_assigned)
table(cls.true)

pairs(x, col = cls.true, pch = 19)
pairs(x, col = labels.salso.VI_assigned, pch = 19)


sim_data.salsoVI <- as.data.frame(cbind(sim_data[, -4], labels.salso.VI_assigned))
colnames(sim_data.salsoVI) <- c("X1", "X2", "X3", "Componente.VI")
sim_data.salsoVI$Componente.VI <- as.factor(sim_data.salsoVI$Componente.VI)

ggplot(sim_data.salsoVI, aes(x = X1, y = X2, color = Componente.VI)) +
  geom_point() +
  stat_ellipse(level = 0.95) + 
  theme_minimal()

plot(summ.VI, type="heatmap")
plot(summ.VI, type="mds") 
plot(summ.VI, type="pairs", data= sim_data[, -3])

loss.salso_VIVI <- partition.loss(sim_data$Componente, labels.salso.VI, loss = VI()) 
loss.salso_VIB <- partition.loss(sim_data$Componente, labels.salso.VI, loss = binder()) 
ari_salso_VI <- adjustedRandIndex(sim_data$Componente, labels.salso.VI) 


save(labels.salso.VI_assigned, loss.salso_VIVI, loss.salso_VIB, ari_salso_VI, time_salso.VI, file= "Greedy_SalsoVI.RData")


# SALSO + loss Binder with a = NULL -------------------------------------------

load("MCMCfit.RData")

library(salso)
library(ggplot2)
library(mclust)


start_time <- Sys.time()
salso.B <- salso(fit$clust, binder(a=NULL)) 
end_time <- Sys.time()

time_salso.B <- difftime(end_time, start_time, units=("secs"))[[1]]


summ.B <- summary(salso.B)

labels.salso.B <- summ.B$estimate

table(labels.salso.B)
table(cls.true)

pairs(x, col = labels.salso.VI, pch = 19)
pairs(x, col = cls.true, pch = 19)


labels.salso.B_assigned <- labels.salso.B

for (j in 1:length(labels.salso.B)){
  if (labels.salso.B[j] == 1) {
    labels.salso.B_assigned[j] <- 3
  } else if (labels.salso.B[j] == 2) {
    labels.salso.B_assigned[j] <- 1
  } else if (labels.salso.B[j] == 3) {
    labels.salso.B_assigned[j] <- 2
  } else if (labels.salso.B[j] == 4) {
    labels.salso.B_assigned[j] <- 5
  } else if (labels.salso.B[j] == 5){
    labels.salso.B_assigned[j] <- 4
  } 
}

table(labels.salso.B_assigned)
table(cls.true)

pairs(x, col = cls.true, pch = 19)
pairs(x, col = labels.salso.B_assigned, pch = 19)


sim_data.salsoB <- as.data.frame(cbind(sim_data[, -4], labels.salso.B_assigned))
colnames(sim_data.salsoB) <- c("X1", "X2", "X3", "Componente.B")
sim_data.salsoB$Componente.B <- as.factor(sim_data.salsoB$Componente.B)

ggplot(sim_data.salsoB, aes(x = X1, y = X2, color = Componente.B)) +
  geom_point() +
  stat_ellipse(level = 0.95) + 
  theme_minimal()


plot(summ.B, type="heatmap")
plot(summ.B, type="mds") 
plot(summ.B, type="pairs", data= sim_data[, -4])

loss.salso_BVI <- partition.loss(sim_data$Componente, labels.salso.B, loss = VI()) 
loss.salso_BB <- partition.loss(sim_data$Componente, labels.salso.B, loss = binder()) 
ari_salso_B <- adjustedRandIndex(sim_data$Componente, labels.salso.B) 


save(labels.salso.B_assigned, loss.salso_BVI, loss.salso_BB, ari_salso_B, time_salso.B, file = "Greedy_SalsoB.RData")


# SALSO + loss Binder con a=list(nClusters=5) -------------------------------------------

load("MCMCfit.RData")

library(salso)
library(ggplot2)
library(mclust)

start_time <- Sys.time()
salso.B2 <- salso(fit$clust, binder(a=list(nClusters=5))) 
end_time <- Sys.time()

time_salso.B2 <- difftime(end_time, start_time, units=("secs"))[[1]]

summ.B2 <- summary(salso.B2)

labels.salso.B2 <- summ.B2$estimate

table(labels.salso.B2)
table(cls.true)

pairs(x, col = labels.salso.B2, pch = 19)
pairs(x, col = cls.true, pch = 19)


labels.salso.B2_assigned <- labels.salso.B2

for (j in 1:length(labels.salso.B2)){
  if (labels.salso.B2[j] == 1) {
    labels.salso.B2_assigned[j] <- 3
  } else if (labels.salso.B2[j] == 2) {
    labels.salso.B2_assigned[j] <- 1
  } else if (labels.salso.B2[j] == 3) {
    labels.salso.B2_assigned[j] <- 2
  } else if (labels.salso.B2[j] == 4) {
    labels.salso.B2_assigned[j] <- 5
  } else {
    labels.salso.B2_assigned[j] <- 4
  }
}

table(labels.salso.B2_assigned)
table(cls.true)

pairs(x, col = cls.true, pch = 19)
pairs(x, col = labels.salso.B2_assigned, pch = 19)


sim_data.salsoB2 <- as.data.frame(cbind(sim_data[, -4], labels.salso.B2_assigned))
colnames(sim_data.salsoB2) <- c("X1", "X2", "X3", "Componente.B2")
sim_data.salsoB2$Componente.B2 <- as.factor(sim_data.salsoB2$Componente.B2)

ggplot(sim_data.salsoB2, aes(x = X1, y = X2, color = Componente.B2)) +
  geom_point() +
  stat_ellipse(level = 0.95) + 
  theme_minimal()

plot(summ.B2, type="heatmap")
plot(summ.B2, type="mds") 
plot(summ.B2, type="pairs", data= sim_data[, -4])

loss.salso_B2VI <- partition.loss(sim_data$Componente, labels.salso.B2, loss = VI()) 
loss.salso_B2B <- partition.loss(sim_data$Componente, labels.salso.B2, loss = binder()) 
ari_salso_B2 <- adjustedRandIndex(sim_data$Componente, labels.salso.B2) 


save(labels.salso.B2_assigned, loss.salso_B2VI, loss.salso_B2B, ari_salso_B2, time_salso.B2, file= "Greedy_SalsoB2.RData")

# Table ----------------------------------------------

load("MCMCfit.RData")
load("Greedy_WadeVI.RData")
load("Greedy_WadeB.RData")
load("Greedy_SalsoVI.RData")
load("Greedy_SalsoB.RData")
load("Greedy_SalsoB2.RData")


library(ggplot2)
library(RColorBrewer)
library(MASS) 
library(reshape)
library(knitr)

kn.wade.VI <- length(unique(labels.wade.minVI_assigned))
kn.wade.B <- length(unique(labels.wade.minbinder_assigned))

kn.salso.VI <-  length(unique(labels.salso.VI_assigned))
kn.salso.B <-  length(unique(labels.salso.B_assigned))
kn.salso.B2 <-  length(unique(labels.salso.B2_assigned))

  
missclass.wade.VI <- n - sum(labels.wade.minVI_assigned == cls.true)
missclass.wade.B <- n - sum(labels.wade.minbinder_assigned == cls.true)
missclass.salso.VI <- n - sum(labels.salso.VI_assigned == cls.true)
missclass.salso.B <- n - sum(labels.salso.B_assigned == cls.true)
missclass.salso.B2 <- n - sum(labels.salso.B2_assigned == cls.true)


loss.wade_VIB
loss.wade_BB
loss.salso_VIB
loss.salso_BB
loss.salso_B2B

loss.wade_VIVI
loss.wade_BVI
loss.salso_VIVI
loss.salso_BVI
loss.salso_B2VI


ari_VI
ari_B
ari_salso_VI
ari_salso_B
ari_salso_B2

Table_perf <- matrix(0,5,6)
rownames(Table_perf) <- c("Greedy_WadeVI", "Greedy_WadeB", "Greedy_SalsoVI", "Greedy_SalsoB", "Greedy_SalsoB2")
colnames(Table_perf) <- c("kn", "Ni", "B(c.hat, ct)", "VI(c.hat, ct)", "ARI", "Tempo Computazionale")



Table_perf[1,c(1:6)] <- round(c(kn.wade.VI, missclass.wade.VI, loss.wade_VIB, loss.wade_VIVI, ari_VI, time_wade.VI), 4)
Table_perf[2,c(1:6)] <- round(c(kn.wade.B, missclass.wade.B, loss.wade_BB, loss.wade_BVI, ari_B, time_wade.B), 4)
Table_perf[3,c(1:6)] <- round(c(kn.salso.VI, missclass.salso.VI, loss.salso_VIB, loss.salso_VIVI, ari_salso_VI, time_salso.VI), 4)
Table_perf[4,c(1:6)] <- round(c(kn.salso.B, missclass.salso.B, loss.salso_BB, loss.salso_BVI, ari_salso_B, time_salso.B), 4)
Table_perf[5,c(1:6)] <- round(c(kn.salso.B2, missclass.salso.B2, loss.salso_B2B, loss.salso_B2VI, ari_salso_B2, time_salso.B2), 4)
Table_perf


# n = 300
#                kn Ni B(c.hat, ct) VI(c.hat, ct)    ARI               Time
# Greedy_WadeVI   5 13       0.0326        0.5323 0.8985              25.0549
# Greedy_WadeB   24 27       0.0384        0.7658 0.8753              32.1067
# Greedy_SalsoVI  5 13       0.0326        0.5323 0.8985               4.2243
# Greedy_SalsoB  16 25       0.0377        0.7310 0.8787             132.0422
# Greedy_SalsoB2  5 14       0.0352        0.5541 0.8909              20.0045


# n = 600
#                kn  Ni B(c.hat, ct) VI(c.hat, ct)    ARI               Time
# Greedy_WadeVI   5  31       0.0391        0.6605 0.8777             117.7035
# Greedy_WadeB   24 279       0.0413        0.8061 0.8666             190.9277
# Greedy_SalsoVI  5  34       0.0428        0.7048 0.8660              16.1122
# Greedy_SalsoB  16  45       0.0419        0.7849 0.8654             341.1415
# Greedy_SalsoB2  5  34       0.0424        0.7217 0.8672              59.5962

# n = 900
#                kn  Ni B(c.hat, ct) VI(c.hat, ct)    ARI               Time
# Greedy_WadeVI   5  48       0.0408        0.7207 0.8724             420.2393
# Greedy_WadeB   32 282       0.0422        0.8815 0.8639             455.5567
# Greedy_SalsoVI  5  52       0.0443        0.7638 0.8614              35.9570
# Greedy_SalsoB  17  67       0.0430        0.8666 0.8622             432.0280
# Greedy_SalsoB2  5  51       0.0433        0.7641 0.8645             151.4777


# n = 1600
#              kn  Ni B(c.hat, ct) VI(c.hat, ct)    ARI                  Time
# Greedy_WadeVI   5  87       0.0420        0.7399 0.8688            3692.0178
# Greedy_WadeB   41 637       0.0408        0.8563 0.8690            2967.4418
# Greedy_SalsoVI  5  84       0.0408        0.7079 0.8725              66.1045
# Greedy_SalsoB  17 110       0.0414        0.8419 0.8678            1194.4606
# Greedy_SalsoB2  5  86       0.0416        0.7328 0.8701             360.4050
 







