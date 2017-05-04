library(Rcpp)
Sys.setenv("PKG_LIBS"=paste(paste0(getwd(), "/Optimise2/cma-es/cma-es/cmaes.o"),
                            paste0(getwd(), "/Optimise2/cma-es/cma-es/boundary_transformation.o")))
sourceCpp("./Optimise2/cma-es/cma-es/cmaes_R.cpp", verbose=TRUE, rebuild=TRUE)
library(ggplot2)

# Simple Optimization Problem


y <- function(x) {
  2 * x * x + 8 * x + 3
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
  stat_function(fun = y) + xlim(-10,10)
p

cost_function <- function(a, b, c) {
  #vals <- seq(from=-10, to=10, length.out = 100)
  #sqrt(sum((y(vals) - (c + b * vals + a * vals * vals))^2))
  
  sum = 0
  for(x in seq(from=-10, to=10, by=0.1)) {
    d = y(x) - a*x*x - b*x - c
    sum <- sum + d*d
  }
  sqrt(sum)
}
cost_function(2, 8, 3)

params <- list()
params$a <- list(min = 1, max = 10) 
params$b <- list(min = 1, max = 10) 
params$c <- list(min = 1, max = 10) 

min <- as.numeric(lapply(params, function(x) x$min))
max <- as.numeric(lapply(params, function(x) x$max))
p <- list(min = min, max = max)
p
#values <- as.numeric(sapply(params, function(x) runif(1, min = x$min, max=x$max)))
values <- c(3, 3, 3)
values
stdDevs <- as.numeric(sapply(params, function(x) (x$max - x$min)/3))
stdDevs


bounds <- boundaryTransformationInit(p)
bounds

set.seed(123)
cmaes <- cmaesInit(values, stdDevs, 123)
cmaes

times <- 100
iterations <- list()
adjustedPop <- list()
for (tm in 1:times) {
  print(tm)
  
  cmaes <- cmaesSamplePopulation(cmaes)
  pop <- cmaes$rgrgx
  
  runs <- cmaes$sp$lambda
  
  #adjustedPop <- list(n = nrow(pop))
  for(l in 1:nrow(pop)) {
    adjustedPop[[l]] <- round(boundaryTransformation(bounds, pop[l,]))
  }

  cost <- vector()
  for (rn in 1:runs) {
    values <- adjustedPop[[rn]]
    
    cost[rn] <- cost_function(round(values[1]), round(values[2]), round(values[3])) 
    
  }
  
  iterations[[tm]] <- cost
  cmaes <- cmaesUpdateDistribution(cmaes, cost)
}  

m <- matrix(unlist(iterations), nrow=times, ncol=nrow(pop), byrow = TRUE)

m <- testConvergence(3)

df <- as.data.frame(m)
names(df) <- c(1:ncol(m))
df$min <- apply(df, 1, min)
df$idx <- as.numeric(rownames(df))

saveRDS(df, "iterations.Rds")
df <- readRDS("iterations.Rds")

p <- ggplot(data = df, aes(x = idx, y = min)) + 
  geom_line() +
  geom_smooth(method = "lm")
p  

library(magrittr)
library(dplyr)
library(tidyr)
df_plot <- df %>%
  select(-min) %>%
  gather(population, cost, -idx)

p <- ggplot() + 
  geom_line(data = df_plot, aes(x = idx, y = cost, group=population, color=population)) +
  xlab('iterations') +
  ylab('Cost')
p

