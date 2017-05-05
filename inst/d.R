setwd("/Users/andreilapanik/Desktop/RocketScience/ltf")
source("testCMA-ES.R")

library(Rcpp)
Sys.setenv("PKG_LIBS"=paste(paste0(getwd(), "/Optimise2/cma-es/cma-es/cmaes.o"),
                            paste0(getwd(), "/Optimise2/cma-es/cma-es/boundary_transformation.o")))
sourceCpp("./Optimise2/cma-es/cma-es/cmaes_R.cpp", verbose=TRUE, rebuild=TRUE)

m <- testConvergence()

source("Tasks_StartMRL_cmaes.r")
