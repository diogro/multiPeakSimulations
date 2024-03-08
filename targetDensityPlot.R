ensure_diag = function(x) (x + t(x))/2
min_dist = 3
space_size = 10
G_obs = ensure_diag(as.matrix(read.csv("data/cov.matrix_Lutreolina.csv", row.names = 1)))

source("./prepareHighDim.R")

df = tidyr::gather(data.frame(beta_enriched = target(n),
                              beta_fit = rbeta(n, shape1 = shapes[[1]][1], 
                                               shape2 = shapes[[1]][2]),
                              random = cor_dist), dist, value)
target_density = ggplot(df, aes(value, group = dist, fill = dist)) + 
  geom_density(alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 1))
save_plot("plots/target_histogram.png", target_density, base_height = 7, base_asp = 1.3)

if(!require(ggridges)){install.packages("ggridges"); library(ggridges)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}
target_density = ggplot(df, aes(value, dist, fill = dist)) + 
  ggridges::geom_density_ridges(alpha = 0.6) + 
  labs(x = expression(paste("Vector correlation between direciton of peak position and ",g[max])),
       y = "Distributions") + 
  scale_fill_viridis(discrete = TRUE) + 
  scale_x_continuous(limits = c(0, 1)) + theme(legend.position = "none") + 
  scale_y_discrete(labels = c("Target\ndistribution", "Fitted Beta", "Observed\nrandom\ndistribution"))
target_density
save_plot("output/plots/target_histogram_ridges.png", target_density, base_height = 7, base_asp = 1.2)

td = target(1000)
hs = hist(td, plot = F, breaks = 5)
library(pander)
intervals = paste(hs$breaks[-6], hs$breaks[-1], sep = "-")
h_table = cbind(c("Correlation\nInterval", "Number\nof vectors"), rbind(intervals,hs$counts))
row.names(h_table) = NULL
pandoc.table(h_table, keep.line.breaks = TRUE)
