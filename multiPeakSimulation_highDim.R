source("./trajectoryTools.R")
source("./plotTools.R")

if(!require(MASS)){install.packages("MASS"); library(MASS)}


library(yamdar)

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(32)

n_peaks = 20
n_traits = 21
space_size = 10


#########################
# G matrices
#########################


# Correlated
G_corr = toadCor
x = eigen(G_corr)$vector[,1]

shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))


peakPool_G_corr_random = randomPeaks(10000, x = eigen(G_corr)$vector[,1], dz_lim = c(3, space_size))
cor_dist = sort(apply(peakPool_G_corr_random, 1, vector_cor, eigen(G_corr)$vector[,1]))



df = tidyr::gather(data.frame(rbeta = rbeta(10000, shape1 = shapes[[1]][1], shape2 = shapes[[1]][2]),
                              vcor = cor_dist), dist, value)
ggplot(df, aes(value, group = dist, fill = dist)) + geom_density(alpha = 0.5)


# Diagonal
G_diag = G_factory(n_traits, rho = 0.1)

peakPool_G_diag_uniform = randomPeaks(10000, x = eigen(G_diag)$vector[,1], intervals = seq(0.05, 1, length.out = 20), prop = rep(1/20, 20), dz_lim = c(3, space_size))
hist(sort(apply(peakPool_G_diag_uniform, 1, vector_cor, eigen(G_diag)$vector[,1])), breaks = 50)

peakPool_G_diag_random = randomPeaks(10000, x = eigen(G_diag)$vector[,1], dz_lim = c(3, space_size))
hist(sort(apply(peakPool_G_diag_random, 1, vector_cor, eigen(G_diag)$vector[,1])), breaks = 50)

peakPool_G_diag_enriched = randomPeaks(10000, x = eigen(G_diag)$vector[,1], intervals = c(0.7, 0.8, 1), prop = c(0.9, 0.06, 0.04), dz_lim = c(3, space_size))
hist(sort(apply(peakPool_G_diag_enriched, 1, vector_cor, eigen(G_diag)$vector[,1])), breaks = 50)


#########################
# Test runs
#########################

runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 4, peakPool = peakPool_G_corr_enriched)
runSimulation("Integrated", G_corr, n_peaks = 20, n_traits, scale = 4, peakPool = peakPool_G_corr_enriched)

#########################
# Simulations
#########################

results_enriched = runTrypitch(G_diag, peakPool_G_corr_enriched, G_corr, peakPool_G_corr_enriched, n = 16, n_peaks = 4)
# results_random   = runTrypitch(G_diag, peakPool_G_corr_random,   G_corr, peakPool_G_corr_random)
# results_uniform  = runTrypitch(G_diag, peakPool_G_corr_uniform,  G_corr, peakPool_G_corr_uniform)
# 
# save(results_enriched,
#      results_random,  
#      results_uniform, file = "plots/results.Rdata")
load("Rdatas/results.Rdata")

plots_enriched = plot_grid(
  plotDzgmax_normdz(results_enriched$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_enriched$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_enriched$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_enriched$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/peakPool_composite_enriched.png", plots_enriched, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)

plots_random = plot_grid(
  plotDzgmax_normdz(results_random$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_random$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_random$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_random$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/peakPool_composite_random.png", plots_enriched, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)

plots_uniform = plot_grid(
  plotDzgmax_normdz(results_uniform$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_uniform$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_uniform$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_uniform$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/peakPool_composite_uniform.png", plots_enriched, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)
