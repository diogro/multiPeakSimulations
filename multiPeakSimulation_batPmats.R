pdf("plots/Rplots.pdf")
source("./plotTools.R")
source("./trajectoryTools.R")

if(!require(MASS)){install.packages("MASS"); library(MASS)}

bat_Ps = list(Ariteus_flavescens = as.matrix(read.csv("./plots/data/cov.matrix_Ariteus_flavescens.csv", row.names = 1)),
              Artibeus_fimbriatus = as.matrix(read.csv("./plots/data/cov.matrix_Artibeus_fimbriatus.csv", row.names = 1)),
              Brachyphylla_cavernarum = as.matrix(read.csv("./plots/data/cov.matrix_Brachyphylla_cavernarum.csv", row.names = 1)))

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(64)

space_size = 10

# G matrices
#########################


## Correlated

# G_corr = G_factory(n_traits, rho = 0.8)
sapply(bat_Ps, CalcR2)
n_traits = dim(bat_Ps[[1]])[1]

n = 100000
peakPool_random = randomPeaks(n, p = n_traits, x = rep(1, n_traits), dz_lim = c(3, space_size))
cor_dist = sort(apply(peakPool_random, 1, vector_cor, eigen(bat_Ps[[1]])$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.15) 

df = tidyr::gather(data.frame(beta_enriched = target(n),
                              beta_fit = rbeta(n, shape1 = shapes[[1]][1], shape2 = shapes[[1]][2]),
                              random = cor_dist), dist, value)
target_density = ggplot(df, aes(value, group = dist, color = dist)) + geom_density(alpha = 0.7) + scale_x_continuous(limits = c(0, 1))
save_plot("plots/target_histogram.png", target_density, base_hei = 7, base_asp = 1.3)

dev.off()
td = target(n)
hs = hist(td, plot = F, breaks = 20)

peakPools = list(randomPeaks(n = 10000, p = n_traits, x = eigen(bat_Ps[[1]])$vector[,1],
                             intervals = hs$breaks[-1], prop = hs$counts/n,dz_lim = c(3, space_size), verbose = FALSE),
                 randomPeaks(n = 10000, p = n_traits, x = eigen(bat_Ps[[2]])$vector[,1],
                             intervals = hs$breaks[-1], prop = hs$counts/n,dz_lim = c(3, space_size), verbose = FALSE),
                 randomPeaks(n = 10000, p = n_traits, x = eigen(bat_Ps[[3]])$vector[,1],
                             intervals = hs$breaks[-1], prop = hs$counts/n,dz_lim = c(3, space_size), verbose = FALSE))

## Diagonal
G_diag = ReplaceDiagonal(G_factory(n_traits, rho = 0.05), diag(G_corr))
eigen(G_diag)$vectors[,1]

cor_dist = sort(apply(peakPool_random, 1, vector_cor, eigen(G_diag)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.25) 
td = target(n)
hs = hist(td, plot = F, breaks = 40)

peakPool_G_diag_enriched = randomPeaks(n = 10000, p = n_traits, x = eigen(G_diag)$vector[,1], intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,dz_lim = c(3, space_size))

# Test runs
#########################

runSimulation("Integrated", bat_Ps[[1]], n_peaks = 1, n_traits, scale = 1, peakPool = peakPools[[1]])
runSimulation("Integrated", bat_Ps[[2]], n_peaks = 1, n_traits, scale = 1, peakPool = peakPools[[2]])
runSimulation("Integrated", bat_Ps[[3]], n_peaks = 1, n_traits, scale = 1, peakPool = peakPools[[3]])

runSimulation("Integrated", bat_Ps[[1]], n_peaks = 50, n_traits, scale = 1, peakPool = peakPools[[1]])
runSimulation("Integrated", bat_Ps[[2]], n_peaks = 50, n_traits, scale = 1, peakPool = peakPools[[2]])
runSimulation("Integrated", bat_Ps[[3]], n_peaks = 50, n_traits, scale = 1, peakPool = peakPools[[3]])

# Simulations
#########################

results_random   = runTrypitch(G_diag, peakPool_random,          G_corr, peakPool_random,          n = 512, n_peaks = 50, scale = 40)
results_enriched = runTrypitchList(bat_Ps, peakPools, n = 32, n_peaks = 10, scale = 1)
# results_uniform  = runTrypitch(G_diag, peakPool_G_corr_uniform,  G_corr, peakPool_G_corr_uniform)
# 
#save(results_enriched,
      #results_random,  
##      results_uniform,
      #file = "plots/results_21d_toadP.Rdata")
# load("Rdatas/results.Rdata")

# Simulations
#########################
                                                          
plots_enriched = plot_grid(
  plotDzgmax_normdz(results_enriched$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_enriched$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_enriched$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_enriched$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/marsupialCovMatrix_composite_enriched.png", plots_enriched, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)

plots_random = plot_grid(
  plotDzgmax_normdz(results_random$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_random$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_random$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_random$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/marsupialCovMatrix_composite_random.png", plots_random, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)

plots_uniform = plot_grid(
  plotDzgmax_normdz(results_uniform$DS, ylim = c(2, space_size), main = "Diagonal G - Single Peak"),
  plotDzgmax_normdz(results_uniform$CS, ylim = c(2, space_size), main = "Integrated G - Single Peak"),
  plotDzgmax_normdz(results_uniform$DM, ylim = c(2, space_size), main = "Diagonal G - Multiple Peaks"),
  plotDzgmax_normdz(results_uniform$CM, ylim = c(2, space_size), main = "Integrated G - Multiple Peaks"),
  ncol = 2, labels = LETTERS[1:4])
save_plot("plots/peakPool_composite_uniform.png", plots_enriched, base_height = 5, base_asp = 1.3, ncol = 2, nrow = 2)
