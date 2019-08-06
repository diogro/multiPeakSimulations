diff_cut_off = 1e-4
max_gens = 10000
max_stand_still = 10
space_size = 6

source("./trajectoryTools.R")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(50)

n_peaks = 20
n_traits = 8

rho = 0.9
while(TRUE){
    G_corr = matrix(rnorm(n_traits*n_traits, rho, 0.05), n_traits, n_traits)
    G_corr = (G_corr + t(G_corr))/2
    diag(G_corr) = rnorm(n_traits, 1, 0.1)
    tryCatch({chol(G_corr); break}, error = function(x) FALSE)
}
chol(G_corr)
eigen(G_corr)
peakPool_G_corr = randomPeaks(10000, x = eigen(G_corr)$vector[,1], steps = 10, min = 2, max = 5)

summary(apply(peakPool_G_corr, 1, Norm))

rho = 0.01
while(TRUE){
    G_diag = matrix(rnorm(n_traits*n_traits, rho, 0.01), n_traits, n_traits)
    G_diag = (G_diag + t(G_diag))/2
    diag(G_diag) = rnorm(n_traits, 1, 0.1)
    tryCatch({chol(G_diag); break}, error = function(x) FALSE)
}
eigen(G_diag)
peakPool_G_diag = randomPeaks(10000, x = eigen(G_diag)$vector[,1], steps = 10, min = 2, max = 5)
summary(apply(peakPool_G_diag, 1, Norm))

# random_peaks = matrix(runif(n_traits*n_peaks, -8, 8), n_peaks, n_traits, byrow = T)
# W_bar_multi = W_bar_factory(random_peaks)
# W_bar_multi_grad = W_bar_gradient_factory(random_peaks)
# W_bar_multi_grad = W_bar_gradient_factory(random_peaks)
# 
# theta_single = matrix(runif(n_traits, -8, 8), 1, n_traits, byrow = T)
# W_bar_single = W_bar_factory(theta_single)
# W_bar_single_grad = W_bar_gradient_factory(theta_single)
# #plotW_bar(W_bar_single)

G = G_corr
p = n_traits
n_peaks = 1
scale = 6
peakPool = peakPool_G_corr
runSimulation = function(G_type = c("Diagonal", "Integrated"), G = NULL,
                         n_peaks = 1, p, rho = 0.7, scale = 6, peakPool){
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
    theta = matrix(peakPool[sample(1:nrow(peakPool), n_peaks),], n_peaks, p)
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
runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 4, peakPool = peakPool_G_corr)
runSimulation("Integrated", G_corr, n_peaks = 50, n_traits, scale = 6, peakPool = peakPool_G_corr)


tic()
results_G_diag_W_single = llply(1:1000, function(x) runSimulation("Diag", G_diag,   1, n_traits, 
                                                                 scale = 2, 
                                                                 peakPool = peakPool_G_diag), 
                                .parallel = TRUE)
results_G_corr_W_single = llply(1:1000, function(x) runSimulation("Inte", G_corr,   1, n_traits, 
                                                                 scale = 4, 
                                                                 peakPool = peakPool_G_corr), 
                                .parallel = TRUE)
results_G_diag_W_multi  = llply(1:1000, function(x) runSimulation("Diag", G_diag, 50, n_traits, 
                                                                 scale = 2, 
                                                                 peakPool = peakPool_G_diag), 
                                .parallel = TRUE)
results_G_corr_W_multi  = llply(1:1000, function(x) runSimulation("Inte", G_corr, 50, n_traits, 
                                                                 scale = 6, 
                                                                 peakPool = peakPool_G_corr), 
                                .parallel = TRUE)
 toc()

save(results_G_diag_W_single,
     results_G_corr_W_single,
     results_G_diag_W_multi,
     results_G_corr_W_multi, file = "results.Rdata")

png("plot_out/peakPool_composite.png", width = 1440, height = 1080)
par(mfrow=c(2, 2))
plotDzgmax_normdz(results_G_diag_W_single, xlim = c(0, 1), ylim = c(1, 7), main = "Diagonal G - Single Peak")
plotDzgmax_normdz(results_G_corr_W_single, xlim = c(0, 1), ylim = c(1, 7), main = "Integrated G - Single Peak")
plotDzgmax_normdz(results_G_diag_W_multi , xlim = c(0, 1), ylim = c(1, 7), main = "Diagonal G - Multiple Peaks")
plotDzgmax_normdz(results_G_corr_W_multi , xlim = c(0, 1), ylim = c(1, 7), main = "Correlated G - Multiple Peaks")
dev.off()


peakPool = 
plot(apply(peakPool[sample(1:10000, 200),], 1, function(x) vector_cor(x, rep(1, 8))))
hist((apply(peakPool_G_diag, 1, Norm)))
prcomp(peakPool)

