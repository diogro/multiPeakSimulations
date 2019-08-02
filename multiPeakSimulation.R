diff_cut_off = 1e-4
max_gens = 10000
max_stand_still = 10
space_size = 10

source("./trajectoryTools.R")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(parallel::detectCores())

n_peaks = 200
n_traits = 8

rho = 0.7
while(TRUE){
    G_corr = matrix(rnorm(n_traits*n_traits, rho, 0.05), n_traits, n_traits)
    G_corr = (G_corr + t(G_corr))/2
    diag(G_corr) = rnorm(n_traits, 1, 0.1)
    tryCatch({chol(G_corr); break}, error = function(x) FALSE)
}
chol(G_corr)
eigen(G_corr)

rho = 0
while(TRUE){
    G_diag = matrix(rnorm(n_traits*n_traits, rho, 0.05), n_traits, n_traits)
    G_diag = (G_diag + t(G_diag))/2
    diag(G_diag) = rnorm(n_traits, 1, 0.1)
    tryCatch({chol(G_diag); break}, error = function(x) FALSE)
}

random_peaks = matrix(runif(n_traits*n_peaks, -8, 8), n_peaks, n_traits, byrow = T)
W_bar_multi = W_bar_factory(random_peaks)
W_bar_multi_grad = W_bar_gradient_factory(random_peaks)
W_bar_multi_grad = W_bar_gradient_factory(random_peaks)

theta_single = matrix(runif(n_traits, -8, 8), 1, n_traits, byrow = T)
W_bar_single = W_bar_factory(theta_single)
W_bar_single_grad = W_bar_gradient_factory(theta_single)
#plotW_bar(W_bar_single)

G = G_corr
p = n_traits
n_peaks = 20
scale = 6
runSimulation = function(G_type = c("Diagonal", "Integrated"), G = NULL,
                         n_peaks = 1, p, rho = 0.7, scale = 6){
    G_type = match.arg(G_type)
    if(is.null(G)){
        if(G_type == "Diagonal"){
            G = diag(p)
            diag(G) = rnorm(p, 1, 0.1)
            gmax = eigen(G)$vectors[,1]
        } else if(G_type == "Integrated"){
            G = matrix(rnorm(p*p, rho, 0.05), p, p)
            G = (G + t(G))/2
            diag(G) = rnorm(p, 1, 0.1)
            gmax = eigen(G)$vectors[,1]
        } else stop("Unknown G type")
    } else
        gmax = eigen(G)$vectors[,1]
    if(n_peaks == 1){
        Surface_type = "Single"
    } else
        Surface_type = "Muliple"
    theta = matrix(runif(p*n_peaks, -space_size, space_size), n_peaks, p, byrow = T)
    W_bar = W_bar_factory(theta)
    W_bar_grad = W_bar_gradient_factory(theta)
    trajectory = calculateTrajectory(rep(0, p), G, W_bar, W_bar_grad, scale = scale)
    trajectory$G_type = G_type
    trajectory$G = G
    trajectory$gmax = gmax
    trajectory$theta = theta
    trajectory$z = trajectory$trajectory[dim(trajectory$trajectory)[1],]
    trajectory$W_bar = W_bar
    trajectory$W_bar_grad = W_bar_grad
    trajectory$Surface_type = Surface_type
    return(trajectory)
}
runSimulation("Integrated", G_corr, n_peaks = 2, n_traits, scale = 10)

results_G_diag_W_single = llply(1:1000, function(x) runSimulation("Diag", G_diag,   1, n_traits, scale = 2), .parallel = TRUE)
results_G_corr_W_single = llply(1:1000, function(x) runSimulation("Inte", G_corr,   1, n_traits, scale = 4), .parallel = TRUE)
results_G_diag_W_multi  = llply(1:1000, function(x) runSimulation("Diag", G_diag, 200, n_traits, scale = 2), .parallel = TRUE)
results_G_corr_W_multi  = llply(1:1000, function(x) runSimulation("Inte", G_corr, 200, n_traits, scale = 6), .parallel = TRUE)

par(mfrow=c(2, 2))
plotDzgmax_normdz(results_G_diag_W_single, xlim = c(0, 1), ylim = c(0, 30), main = "Diagonal G - Single Peak")
plotDzgmax_normdz(results_G_corr_W_single, xlim = c(0, 1), ylim = c(0, 30), main = "Integrated G - Single Peak")
plotDzgmax_normdz(results_G_diag_W_multi , xlim = c(0, 1), ylim = c(0, 30), main = "Diagonal G - Multiple Peaks")
plotDzgmax_normdz(results_G_corr_W_multi , xlim = c(0, 1), ylim = c(0, 30), main = "Correlated G - Multiple Peaks")
