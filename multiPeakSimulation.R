diff_cut_off = 1e-4
max_gens = 10000
max_stand_still = 10
space_size = 10

source("./trajectoryTools.R")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(parallel::detectCores())

n_peaks = 100
n_traits = 6

rho = 0.8
G_corr = matrix(rnorm(n_traits*n_traits, rho, 0.05), n_traits, n_traits)
G_corr = (G_corr + t(G_corr))/2
diag(G_corr) = rnorm(n_traits, 1, 0.1)
gmax_corr = eigen(G_corr)$vectors[,1]

G_diag = diag(n_traits)
diag(G_diag) = rnorm(n_traits, 1, 0.1)
gmax_diag = eigen(G_diag)$vectors[,1]

random_peaks = matrix(runif(n_traits*n_peaks, -8, 8), n_peaks, n_traits, byrow = T)
W_bar_multi = W_bar_factory(random_peaks)
W_bar_multi_grad = W_bar_gradient_factory(random_peaks)
W_bar_multi_grad = W_bar_gradient_factory(random_peaks)

theta_single = matrix(runif(n_traits, -8, 8), 1, n_traits, byrow = T)
W_bar_single = W_bar_factory(theta_single)
W_bar_single_grad = W_bar_gradient_factory(theta_single)
#plotW_bar(W_bar_single)


n_sims = 1000
random_start = matrix(runif(n_traits*n_sims, -space_size, space_size),
                      n_sims, n_traits, byrow = T)
results_G_diag_W_single = alply(random_start, 1, calculateTrajectory, G_diag, W_bar_single, W_bar_single_grad, scale = 6, .progress = "text")
results_G_corr_W_single = alply(random_start, 1, calculateTrajectory, G_corr, W_bar_single, W_bar_single_grad, scale = 6, .progress = "text")
results_G_diag_W_multi  = alply(random_start, 1, calculateTrajectory, G_diag, W_bar_multi,  W_bar_multi_grad,  scale = 6, .progress = "text")
results_G_corr_W_multi  = alply(random_start, 1, calculateTrajectory, G_corr, W_bar_multi,  W_bar_multi_grad,  scale = 6, .progress = "text")


plot(laply(results_G_diag_W_single, function(x) x$trajectory[dim(x$trajectory)[1],]))
plot(laply(results_G_corr_W_single, function(x) x$trajectory[dim(x$trajectory)[1],]))
plot(laply(results_G_diag_W_multi , function(x) x$trajectory[dim(x$trajectory)[1],]))
plot(laply(results_G_corr_W_multi , function(x) x$trajectory[dim(x$trajectory)[1],]))

plotDzgmax_normdz(results_G_diag_W_single, gmax_diag)
plotDzgmax_normdz(results_G_corr_W_single, gmax_corr)
plotDzgmax_normdz(results_G_diag_W_multi , gmax_diag)
plotDzgmax_normdz(results_G_corr_W_multi , gmax_corr)
