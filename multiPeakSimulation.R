
source("./trajectoryTools.R")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(32)

n_peaks = 20
n_traits = 8
space_size = 6


#########################
# G matrices
#########################


# Correlated
G_corr = G_factory(n_traits, rho = 0.9)

peakPool_G_corr_uniform = randomPeaks(10000, x = eigen(G_corr)$vector[,1], intervals = seq(0.025, 1, length.out = 40), prop = rep(1/40, 40), dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_corr_uniform, 1, vector_cor, eigen(G_corr)$vector[,1])), breaks = 50)

peakPool_G_corr_random = randomPeaks(10000, x = eigen(G_corr)$vector[,1], dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_corr_random, 1, vector_cor, eigen(G_corr)$vector[,1])), breaks = 50)

peakPool_G_corr_enriched = randomPeaks(10000, x = eigen(G_corr)$vector[,1], intervals = c(0.7, 0.8, 1), prop = c(0.9, 0.06, 0.04), dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_corr_enriched, 1, vector_cor, eigen(G_corr)$vector[,1])), breaks = 50)


# Diagonal
G_diag = G_factory(n_traits, rho = 0.1)

peakPool_G_diag_uniform = randomPeaks(10000, x = eigen(G_diag)$vector[,1], intervals = seq(0.025, 1, length.out = 40), prop = rep(1/40, 40), dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_diag_uniform, 1, vector_cor, eigen(G_diag)$vector[,1])), breaks = 50)

peakPool_G_diag_random = randomPeaks(10000, x = eigen(G_diag)$vector[,1], dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_diag_random, 1, vector_cor, eigen(G_diag)$vector[,1])), breaks = 50)

peakPool_G_diag_enriched = randomPeaks(10000, x = eigen(G_diag)$vector[,1], intervals = c(0.7, 0.8, 1), prop = c(0.9, 0.06, 0.04), dz_lim = c(2, 6))
hist(sort(apply(peakPool_G_diag_enriched, 1, vector_cor, eigen(G_diag)$vector[,1])), breaks = 50)


#########################
# Test runs
#########################

runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 4, peakPool = peakPool_G_corr)
runSimulation("Integrated", G_corr, n_peaks = 20, n_traits, scale = 4, peakPool = peakPool_G_corr)

#########################
# Simulations
#########################

runTrypitch = function(G_diag, peakPool_diag, G_corr, peakPool_corr, n = 1000, n_peaks = 50, parallel = TRUE){
  results_G_diag_W_single = llply(1:n, function(x) runSimulation("Diag", G_diag,   1, n_traits, 
                                                                    scale = 4, 
                                                                    peakPool = peakPool_diag), 
                                  .parallel = parallel)
  results_G_corr_W_single = llply(1:n, function(x) runSimulation("Inte", G_corr,   1, n_traits, 
                                                                    scale = 4, 
                                                                    peakPool = peakPool_corr), 
                                  .parallel = parallel)
  results_G_diag_W_multi  = llply(1:n, function(x) runSimulation("Diag", G_diag, n_peaks, n_traits, 
                                                                    scale = 4, 
                                                                    peakPool = peakPool_diag), 
                                  .parallel = parallel)
  results_G_corr_W_multi  = llply(1:n, function(x) runSimulation("Inte", G_corr, n_peaks, n_traits, 
                                                                    scale = 4, 
                                                                    peakPool = peakPool_corr), 
                                  .parallel = parallel)
  list(DS = results_G_diag_W_single,
       CS = results_G_corr_W_single,
       DM = results_G_diag_W_multi,
       CM = results_G_corr_W_multi)
}

results_enriched = runTrypitch(G_diag, peakPool_G_diag_enriched, G_corr, peakPool_G_corr_enriched, n=2, n_peaks = 5, parallel = FALSE)
results_random   = runTrypitch(G_diag, peakPool_G_diag_random,   G_corr, peakPool_G_corr_random)
results_uniform  = runTrypitch(G_diag, peakPool_G_diag_uniform,  G_corr, peakPool_G_corr_uniform)


#save(results_G_diag_W_single,
     #results_G_corr_W_single,
     #results_G_diag_W_multi,
     #results_G_corr_W_multi, file = "results.Rdata")

png("plot_out/peakPool_composite_enriched.png", width = 1440, height = 1080)
par(mfrow=c(2, 2))
plotDzgmax_normdz(results_enriched$DS, xlim = c(0, 1), ylim = c(1, 6), main = "Diagonal G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_enriched$CS, xlim = c(0, 1), ylim = c(1, 6), main = "Integrated G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_enriched$DM, xlim = c(0, 1), ylim = c(1, 6), main = "Diagonal G - Multiple Peaks")
abline(h=2)
plotDzgmax_normdz(results_enriched$CM, xlim = c(0, 1), ylim = c(1, 6), main = "Integrated G - Multiple Peaks")
abline(h=2)
dev.off()

png("plot_out/peakPool_composite_random.png", width = 1440, height = 1080)
par(mfrow=c(2, 2))
plotDzgmax_normdz(results_random$DS, xlim = c(0, 1), ylim = c(1, 6), main = "Diagonal G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_random$CS, xlim = c(0, 1), ylim = c(1, 6), main = "Integrated G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_random$DM, xlim = c(0, 1), ylim = c(1, 6), main = "Diagonal G - Multiple Peaks")
abline(h=2)
plotDzgmax_normdz(results_random$CM, xlim = c(0, 1), ylim = c(1, 6), main = "Integrated G - Multiple Peaks")
abline(h=2)
dev.off()

png("plot_out/peakPool_composite_uniform.png", width = 1440, height = 1080)
par(mfrow=c(2, 2))
plotDzgmax_normdz(results_uniform$DS, xlim = c(0, 1), ylim = c(1, 6), main = "Diagonal G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_uniform$CS, xlim = c(0, 1), ylim = c(1, 6), main = "Integrated G - Single Peak")
abline(h=2)
plotDzgmax_normdz(results_uniform$DM, xlim = c(0, 1), ylim = c(1, 6), main = "Diagonal G - Multiple Peaks")
abline(h=2)
plotDzgmax_normdz(results_uniform$CM, xlim = c(0, 1), ylim = c(1, 6), main = "Integrated G - Multiple Peaks")
abline(h=2)
dev.off()
