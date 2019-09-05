source("./trajectoryTools.R")
source("./plotTools.R")

if(!require(MASS)){install.packages("MASS"); library(MASS)}

devtools::install_github("diogro/yamda-r", subdir = "package")
library(yamdar)

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(32)

n_peaks = 20
n_traits = 40 
space_size = 10


#########################
# G matrices
#########################

rbeta_mixture = function(n, shapes1, shapes2, alpha){
  f1 = function(x) rbeta(x, shapes1[1], shapes1[2])
  f2 = function(x) rbeta(x, shapes2[1], shapes2[2])
  out = numeric(n)
  for(i in 1:n){
    if(runif(1) > alpha){
      out[i] = f1(1)
    }else 
      out[i] = f2(1)
  }
  out
}

# Correlated
G_corr = G_factory(n_traits, rho = 0.1)
G_corr = toadCor
n_traits = dim(toadCor)[1]
x = eigen(G_corr)$vector[,1]

peakPool_G_corr_random = randomPeaks(10000, p = n_traits, x = eigen(G_corr)$vector[,1], dz_lim = c(3, space_size))
cor_dist = sort(apply(peakPool_G_corr_random, 1, vector_cor, eigen(G_corr)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

df = tidyr::gather(data.frame(rbeta_En = rbeta_mixture(10000, shapes[[1]], c(2, 3), 0.15),
                              rbeta = rbeta(10000, shape1 = shapes[[1]][1], shape2 = shapes[[1]][2]),
                              vcor = cor_dist), dist, value)
ggplot(df, aes(value, group = dist, color = dist)) + geom_density(alpha = 0.7)

td = rbeta_mixture(n, shapes[[1]], c(2, 3), 0.15)
hs = hist(td, plot = T, breaks = 40)

peakPool_G_corr_enriched = randomPeaks(n = 10000, p = n_traits, x = eigen(G_corr)$vector[,1], intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,dz_lim = c(3, space_size))
cor_dist_2 = sort(apply(peakPool_G_corr_enriched, 1, vector_cor, eigen(G_corr)$vector[,1]))
hist(cor_dist_2)

# Diagonal
G_diag = G_factory(n_traits, rho = 0.1)

peakPool_G_diag_random = randomPeaks(10000, p = n_traits, x = eigen(G_diag)$vector[,1], dz_lim = c(3, space_size))
cor_dist = sort(apply(peakPool_G_diag_random, 1, vector_cor, eigen(G_diag)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

df = tidyr::gather(data.frame(rbeta_En = rbeta_mixture(10000, shapes[[1]], c(2, 3), 0.15),
                              rbeta = rbeta(10000, shape1 = shapes[[1]][1], shape2 = shapes[[1]][2]),
                              vcor = cor_dist), dist, value)
ggplot(df, aes(value, group = dist, color = dist)) + geom_density(alpha = 0.7)

td = rbeta_mixture(n, shapes[[1]], c(2, 3), 0.15)
hs = hist(td, plot = T, breaks = 40)

peakPool_G_diag_enriched = randomPeaks(n = 10000, p = n_traits, x = eigen(G_diag)$vector[,1], intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,dz_lim = c(3, space_size))


#########################
# Test runs
#########################

runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 10, peakPool = peakPool_G_corr_enriched)
runSimulation("Integrated", G_corr, n_peaks = 20, n_traits, scale = 10, peakPool = peakPool_G_corr_enriched)

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
