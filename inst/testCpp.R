library(Rcpp)
Sys.setenv("PKG_LIBS"=paste(paste0(getwd(), "/Optimise2/cma-es/cma-es/cmaes.o"),
                            paste0(getwd(), "/Optimise2/cma-es/cma-es/boundary_transformation.o")))
sourceCpp("./Optimise2/cma-es/cma-es/cmaes_R.cpp", verbose=TRUE, rebuild=TRUE)
sp = list(filename = "test", 
         flgsupplemented = 1,
         N = 3, 
         seed = 2,
         xstart=c(1.5, 2.5, 3.5),
         typicalX=c(2.2, 3.3, 4.4)
         )
testCmaesReadPara(sp, 3)

rand_in = list(startseed = 1,
            aktseed = 2,
            aktrand = 3,
            rgrand = c(1:32),
            flgstored = 4,
            hold = 1.234)

rand_out = testCmaesRand(rand_in)
rand_out

mat_in = matrix(c(1:9), nrow=3, ncol=3)
mat_out = testMatrix(mat_in, 3, 3) 
mat_out

params <- vector(length = 11)
params[1] <- round(runif(1, min=1    , max=50   )) #VOLATILITY_MULTIPLIER                                                                                #####
params[2] <- round(runif(1, min=1    , max=12   )) #VOLATILITY_PERIODS                                                                               #####
params[3] <- round(runif(1, min=1    , max=12   )) #AVERAGE_PERIODS                                                                                #####
params[4] <- runif(1, min=0.3  , max=0.9  ) #SMOOTHING_COEFFICIENT                                                                                 #####
params[5] <- round(runif(1, min=1    , max=10   ))  #PAST_DATAPOINTS                                                                              #####
params[6] <- runif(1, min=0.001, max=0.01 ) #NOISE_THRESHOLD                                                                                #####
params[7] <- round(runif(1, min=2    , max=10   )) #CONSECUTIVE_CLOSES                                                                               #####
params[8] <- runif(1, min=0    , max=0.01 ) #LOSS_TARGET                                                                                #####
params[9] <- runif(1, min=0.001, max=0.008) #PROFIT_TARGET                                                                                #####
params[10] <- round(runif(1, min=2    , max=30   )) #HOLDING_PERIOD                                                                                #####
params[11] <- round(runif(1, min=1    , max=100  ))  #SLEEPING_PERIOD 

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

sourceCpp("./Optimise2/cma-es/cma-es/cmaes_R.cpp", verbose=TRUE, rebuild=TRUE)

bounds <- boundaryTransformationInit(p)
bounds

cmaes <- cmaesInit(values, stdDevs)
cmaes

cmaes1 <- cmaesInitWithSamplePopulation(values, stdDevs)
cmaes1$rgrgx

cmaes_out <- testCmaes(cmaes)
all.equal(cmaes_out, cmaes)

cmaes2 <- cmaesSamplePopulation(cmaes)
cmaes2

all.equal(cmaes1, cmaes1)

values <- cmaes$rgrgx[1,]
values

pop <- boundaryTransformation(bounds, values)
pop

rgFunVal <- c(1:6)

cmaes <- cmaesUpdateDistribution(cmaes, rgFunVal)
cmaes

