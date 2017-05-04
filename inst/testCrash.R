library(Rcpp)
Sys.setenv("PKG_LIBS"=paste(paste0(getwd(), "/Optimise2/cma-es/cma-es/cmaes.o"),
                            paste0(getwd(), "/Optimise2/cma-es/cma-es/boundary_transformation.o")))
cmaes_R <- paste0(getwd(), "/Optimise2/cma-es/cma-es/cmaes_R.cpp")
sourceCpp(cmaes_R, verbose=TRUE, rebuild=TRUE)

params <- vector("list", 2)

params[[1]]$name <- "VOLATILITY_MULTIPLIER"
params[[1]]$value <- round(runif(1, min=1, max=50))
params[[1]]$min <- 1
params[[1]]$max <- 50

params[[2]]$name <- "VOLATILITY_PERIODS"
params[[2]]$value <- round(runif(1, min=1, max=12))
params[[2]]$min <- 1
params[[2]]$max <- 12

min <- as.numeric(lapply(params, function(x) x$min))
max <- as.numeric(lapply(params, function(x) x$max))
p <- list(min = min, max = max)
p

values <- as.numeric(sapply(params, function(x) x$value))
values
stdDevs <- as.numeric(sapply(params, function(x) (x$max - x$min)/3))
stdDevs



cmaes <- cmaesInit(values, stdDevs)
cmaes

pop <- cmaesSamplePopulation(cmaes)
pop
