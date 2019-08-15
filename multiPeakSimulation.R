
source("./trajectoryTools.R")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(32)

n_peaks = 20
n_traits = 8
space_size = 6


#########################
# G matrices
#########################


G_corr = G_factory(n_traits, rho = 0.9)

peakPool_G_corr = randomPeaks(10000, x = eigen(G_corr)$vector[,1], intervals = c(0.7, 0.8, 1), prop = c(0.9, 0.06, 0.04), dz_lim = c(2, 6))
#peakPool_G_corr = randomPeaks(10000, x = eigen(G_corr)$vector[,1], dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_corr, 1, vector_cor, eigen(G_corr)$vector[,1])), breaks = 50)

G_diag = G_factory(n_traits, rho = 0.1)

peakPool_G_diag = randomPeaks(10000, x = eigen(G_diag)$vector[,1], intervals = c(0.7, 0.8, 1), prop = c(0.9, 0.06, 0.04), dz_lim = c(2, 6))
#peakPool_G_diag = randomPeaks(10000, x = eigen(G_diag)$vector[,1], dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_diag, 1, vector_cor, eigen(G_diag)$vector[,1])), breaks = 50)

#########################
# Test runs
#########################

runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 4, peakPool = peakPool_G_corr)
runSimulation("Integrated", G_corr, n_peaks = 20, n_traits, scale = 4, peakPool = peakPool_G_corr)

#########################
# Simulations
#########################

tic()
results_G_diag_W_single = llply(1:1000, function(x) runSimulation("Diag", G_diag,   1, n_traits, 
                                                                 scale = 4, 
                                                                 peakPool = peakPool_G_diag), 
                                .parallel = TRUE)
results_G_corr_W_single = llply(1:1000, function(x) runSimulation("Inte", G_corr,   1, n_traits, 
                                                                 scale = 4, 
                                                                 peakPool = peakPool_G_corr), 
                                .parallel = TRUE)
results_G_diag_W_multi  = llply(1:1000, function(x) runSimulation("Diag", G_diag, 6, n_traits, 
                                                                 scale = 4, 
                                                                 peakPool = peakPool_G_diag), 
                                .parallel = TRUE)
results_G_corr_W_multi  = llply(1:1000, function(x) runSimulation("Inte", G_corr, 50, n_traits, 
                                                                 scale = 4, 
                                                                 peakPool = peakPool_G_corr), 
                                .parallel = TRUE)
toc()

#save(results_G_diag_W_single,
     #results_G_corr_W_single,
     #results_G_diag_W_multi,
     #results_G_corr_W_multi, file = "results.Rdata")

png("plot_out/peakPool_composite.png", width = 1440, height = 1080)
par(mfrow=c(2, 2))
plotDzgmax_normdz(results_G_diag_W_single, xlim = c(0, 1), ylim = c(2, 6), main = "Diagonal G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_G_corr_W_single, xlim = c(0, 1), ylim = c(2, 6), main = "Integrated G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_G_diag_W_multi , xlim = c(0, 1), ylim = c(1, 6), main = "Diagonal G - Multiple Peaks")
abline(h=2)
plotDzgmax_normdz(results_G_corr_W_multi , xlim = c(0, 1), ylim = c(1, 6), main = "Integrated G - Multiple Peaks")
abline(h=2)
dev.off()

