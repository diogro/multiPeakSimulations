pdf("plots/Rplots.pdf")
source("./plotTools.R")
source("./trajectoryTools.R")

if(!require(MASS)){install.packages("MASS"); library(MASS)}

load("./orders.Rdata")

if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(4)

space_size = 10

# G matrices
#########################


## Correlated

# G_corr = G_factory(n_traits, rho = 0.8)
mammal.orders = llply(mammal.orders, function(x) (x + t(x))/2)
max_int = max(sapply(mammal.orders, CalcEigenVar))
G_obs = mammal.orders$Lutreolina
n_traits = dim(G_obs)[1]

set.seed(42)
Random_G = RandomMatrix(n_traits, variance = diag(G_obs), LKJ = FALSE)
G_corr = expEigenVal(Random_G, 4)
CalcEigenVar(G_corr)
x = eigen(G_corr)$vector[,1]

n = 100000
peakPool_random = randomPeaks(n, p = n_traits, 
                              x = eigen(G_corr)$vector[,1], 
                              dz_lim = c(3, space_size))

cor_dist = sort(apply(peakPool_random, 1, vector_cor, eigen(G_corr)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))
target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.15) 

df = tidyr::gather(data.frame(beta_enriched = target(n),
                              beta_fit = rbeta(n, shape1 = shapes[[1]][1], 
                                               shape2 = shapes[[1]][2]),
                              random = cor_dist), dist, value)
target_density = ggplot(df, aes(value, group = dist, fill = dist)) + 
  geom_density(alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 1))
save_plot("plots/target_histogram.png", target_density, 
          base_height = 7, base_asp = 1.3)

if(!require(ggridges)){install.packages("ggridges"); library(ggridges)}
target_density = ggplot(df, aes(value, dist, fill = dist)) + 
  ggridges::geom_density_ridges(alpha = 0.6) + 
  scale_x_continuous(limits = c(0, 1))
target_density

plot(diag(G_corr))
dev.off()
td = target(n)
hs = hist(td, plot = F, breaks = 20)

peakPool_G_corr_enriched = randomPeaks(n = 10000, 
                                       p = n_traits, 
                                       x = eigen(G_corr)$vector[,1], 
                                       intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,
                                       dz_lim = c(3, space_size), verbose = FALSE)

## Diagonal

G_diag = expEigenVal(Random_G, 0.2)
CalcEigenVar(G_diag)
eigen(G_diag)$vectors[,1]

cor_dist = sort(apply(peakPool_random, 1, vector_cor, eigen(G_diag)$vector[,1]))
shapes = fitdistr(cor_dist, dbeta, list(shape1=1, shape2=10))

target = function(x) rbeta_mixture(x, shapes[[1]], c(1, 1), 0.15) 
td = target(n)
hs = hist(td, plot = F, breaks = 20)

peakPool_G_diag_enriched = randomPeaks(n = 10000, 
                                       p = n_traits, 
                                       x = eigen(G_diag)$vector[,1], 
                                       intervals = hs$breaks[-1], 
                                       prop = hs$counts/n,
                                       dz_lim = c(3, space_size))

# Test runs
#########################

runSimulation("Integrated", G_corr, n_peaks = 1, n_traits, scale = 40, peakPool = peakPool_G_corr_enriched)
runSimulation("Diagonal", G_diag, n_peaks = 1, n_traits, scale = 40, peakPool = peakPool_G_diag_enriched)
runSimulation("Integrated", G_corr, n_peaks = 50, n_traits, scale = 40, peakPool = peakPool_random)
runSimulation("Diagonal", G_diag, n_peaks = 10, n_traits, scale = 40, peakPool = peakPool_random)

# Simulations
#########################

results_random   = runTrypitch(G_diag, peakPool_random,          
                               G_corr, peakPool_random,          
                               n = 256, n_peaks = 50, scale = 40)
results_enriched = runTrypitch(G_diag, peakPool_G_diag_enriched, 
                               G_corr, peakPool_G_corr_enriched, 
                               n = 256, n_peaks = 50, scale = 40)
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
