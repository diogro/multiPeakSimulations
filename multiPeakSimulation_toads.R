pdf("plots/Rplots.pdf")
source("./plotTools.R")
source("./trajectoryTools.R")

if(!require(MASS)){install.packages("MASS"); library(MASS)}

#devtools::install_github("diogro/yamda-r", subdir = "package")
library(yamdar)

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(64)

space_size = 10

# G matrices
#########################


## Correlated

# G_corr = G_factory(n_traits, rho = 0.8)
G_corr = toadCor
n_traits = dim(toadCor)[1]
x = eigen(G_corr)$vector[,1]

n = 100000
peakPool_G_corr_random = randomPeaks(n, p = n_traits, x = eigen(G_corr)$vector[,1], dz_lim = c(3, space_size))
cor_dist = sort(apply(peakPool_G_corr_random, 1, vector_cor, eigen(G_corr)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.25) 

df = tidyr::gather(data.frame(beta_enriched = target(n),
                              beta_fit = rbeta(n, shape1 = shapes[[1]][1], shape2 = shapes[[1]][2]),
                              random = cor_dist), dist, value)
target_density = ggplot(df, aes(value, group = dist, color = dist)) + geom_density(alpha = 0.7) + scale_x_continuous(limits = c(0, 1))
save_plot("plots/target_histogram.png", target_density, base_hei = 7, base_asp = 1.3)

td = target(n)
hs = hist(td, plot = F, breaks = 40)

peakPool_G_corr_enriched = randomPeaks(n = 10000, p = n_traits, x = eigen(G_corr)$vector[,1], intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,dz_lim = c(3, space_size))

## Diagonal
G_diag = G_factory(n_traits, rho = 0.05)
eigen(G_diag)$vectors[,1]

peakPool_G_diag_random = randomPeaks(10000, p = n_traits, x = eigen(G_diag)$vector[,1], dz_lim = c(3, space_size))
cor_dist = sort(apply(peakPool_G_diag_random, 1, vector_cor, eigen(G_diag)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.25) 
td = target(n)
hs = hist(td, plot = F, breaks = 40)

peakPool_G_diag_enriched = randomPeaks(n = 10000, p = n_traits, x = eigen(G_diag)$vector[,1], intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,dz_lim = c(3, space_size))

# Test runs
#########################

runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 10, peakPool = peakPool_G_corr_enriched)
runSimulation("Diagonal", G_diag, n_peaks = 1, n_traits, scale = 10, peakPool = peakPool_G_corr_enriched)
runSimulation("Integrated", G_corr, n_peaks = 10, n_traits, scale = 10, peakPool = peakPool_G_corr_enriched)
runSimulation("Diagonal", G_diag, n_peaks = 10, n_traits, scale = 10, peakPool = peakPool_G_corr_enriched)

# Simulations
#########################
results_enriched = runTrypitch(G_diag, peakPool_G_corr_enriched, G_corr, peakPool_G_corr_enriched, n = 1024, n_peaks = 50, scale = 10)
results_random   = runTrypitch(G_diag, peakPool_G_corr_random,   G_corr, peakPool_G_corr_random, n = 1024, n_peaks = 50, scale = 10)
# results_uniform  = runTrypitch(G_diag, peakPool_G_corr_uniform,  G_corr, peakPool_G_corr_uniform)
# 
 save(results_enriched,
      results_random,  
#      results_uniform,
      file = "plots/results_21d_toadP.Rdata")
# load("Rdatas/results.Rdata")

# Simulations
#########################
                                                          
plots_enriched = plot_grid(
  plotDzgmax_normdz(results_enriched$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_enriched$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_enriched$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_enriched$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/toadMatrix_composite_enriched.png", plots_enriched, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)

plots_random = plot_grid(
  plotDzgmax_normdz(results_random$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_random$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_random$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_random$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/toadMatrix_composite_random.png", plots_random, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)

plots_uniform = plot_grid(
  plotDzgmax_normdz(results_uniform$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_uniform$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_uniform$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_uniform$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/peakPool_composite_uniform.png", plots_enriched, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)
