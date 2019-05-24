diff_cut_off = 1e-1
max_gens = 1000
max_stand_still = 3
space_size = 10

source("./trajectoryTools.R")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(parallel::detectCores())

n_peaks = 1000
n_traits = 5
random_peaks = matrix(runif(n_traits*n_peaks, -8, 8), n_peaks, n_traits, byrow = T)
mean(dist(random_peaks))
W_bar_multi = W_bar_factory(random_peaks)
#plotW_bar(W_bar_multi)

theta_single = matrix(c(3, 3), 1, 2, byrow = T)
W_bar_single = W_bar_factory(theta_single)
plotW_bar(W_bar_single)

rho = 0.8
G_corr = matrix(rho, n_traits, n_traits)
diag(G_corr) = rnorm(n_traits, 1, 0.1)
gmax_corr = eigen(G_corr)$vectors[,1]

G_diag = diag(n_traits)
diag(G_diag) = rnorm(n_traits, 1, 0.1)
gmax_diag = eigen(G_diag)$vectors[,1]


n_sims = 100
random_start = matrix(runif(n_traits*n_sims, -space_size, space_size),
                      n_sims, n_traits, byrow = T)
results_corr = alply(random_start, 1, calculateTrajectory, G_corr, W_bar_multi, scale = 6, .parallel=TRUE)


plot(laply(results_corr, function(x) x$trajectory[dim(x$trajectory)[1],]))
plot(laply(results_corr, function(x) x$net_dz))

plotDzgmax_normdz(results_corr, gmax_corr)
plotDzgmax_normdz(results_diag, gmax_diag)
