source("./plotTools.R")
source("./trajectoryTools.R")

n_traits = dim(G_obs)[1]

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
if(diag_mat_type == "low"){
  G_diag = expEigenVal(G_obs, 0.2)
} else if(diag_mat_type == "diag"){
  G_diag = diag(diag(G_obs))
} else{
  log4r_error("Unknown low integration matrix type.")
  stop()
}
log4r_info(paste0("Eigen variance in observed matrix: ", CalcEigenVar(G_obs)))
log4r_info(paste0("Eigen variance in unconstrained matrix: ", CalcEigenVar(G_diag)))

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
