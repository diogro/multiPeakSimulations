source("./trajectoryTools.R")

w_cov = matrix(c(1.0, 0.0,
                 0.0, 1.0), ncol = 2)


W_bar_multi = function(x) {
  log(
    dmvnorm(x, mean = c(4, 5), sigma = w_cov) +
    dmvnorm(x, mean = c(1, 5), sigma = w_cov) +
    dmvnorm(x, mean = c(2, 2), sigma = w_cov) +
    dmvnorm(x, mean = c(7, 5), sigma = w_cov) +
    dmvnorm(x, mean = c(5, 2), sigma = w_cov))
}
plotW_bar(W_bar_multi)

W_bar_single = function(x) {
  log(dmvnorm(x, mean = c(3, 3), sigma = w_cov))
}

G_corr = matrix(c(1.1, 0.8,
             0.8, 1.0), ncol = 2)
gmax_corr = eigen(G_corr)$vectors[,1]

G_diag = matrix(c(1.1, 0.0,
                  0.0, 1.0), ncol = 2)
gmax_diag = eigen(G_diag)$vectors[,1]


n_sims = 10
results_corr = vector("list", n_sims)
for(i in 1:n_sims){
  random_start = runif(2, -space_size, space_size)
  results_corr[[i]] = calculateTrajectory(G_corr, W_bar_multi, start_position = random_start, scale = 2)
}

plotW_bar(W_bar_multi)
for(i in 1:n_sims){
  points(results_corr[[i]]$trajectory-2, col=mypalette[i])  
}


plotDzgmax_normdz(results_corr, gmax_corr)
plotDzgmax_normdz(results_diag, gmax_diag)
