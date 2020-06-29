source("./plotTools.R")
source("./trajectoryTools.R")

load("./orders.Rdata")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
n_cores = min(detectCores()-1, 64)
registerDoMC(n_cores)

min_dist = 3
space_size = 10

# G matrices
#########################

## Correlated

# G_corr = G_factory(n_traits, rho = 0.8)
mammal.orders = llply(mammal.orders, function(x) (x + t(x))/2)
max_int = max(sapply(mammal.orders, CalcEigenVar))
G_obs = mammal.orders$Lutreolina
n_traits = dim(G_obs)[1]

set.seed(42)
#Random_G = RandomMatrix(n_traits, variance = diag(G_obs), LKJ = FALSE)
#G_corr = expEigenVal(Random_G, 4)
G_corr = G_obs
CalcEigenVar(G_corr)
x = eigen(G_corr)$vector[,1]

n = 100000
peakPool_random = randomPeaks(n, p = n_traits, 
                              x = eigen(G_corr)$vector[,1], 
                              dz_lim = c(min_dist, space_size))

cor_dist = sort(apply(peakPool_random, 1, vector_cor, eigen(G_corr)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))
target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.15) 

td = target(n)
hs = hist(td, plot = F, breaks = 30)

peakPool_G_corr_enriched = randomPeaks(n = 10000, 
                                       p = n_traits, 
                                       x = eigen(G_corr)$vector[,1], 
                                       intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,
                                       dz_lim = c(min_dist, space_size), verbose = FALSE)


## Diagonal
G_diag = expEigenVal(G_obs, 0.2)
CalcEigenVar(G_diag)
eigen(G_diag)$vectors[,1]

cor_dist = sort(apply(peakPool_random, 1, vector_cor, eigen(G_diag)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.15) 
td = target(n)
hs = hist(td, plot = F, breaks = 20)

peakPool_G_diag_enriched = randomPeaks(n = 10000, 
                                       p = n_traits, 
                                       x = eigen(G_diag)$vector[,1], 
                                       intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,
                                       dz_lim = c(min_dist, space_size))
