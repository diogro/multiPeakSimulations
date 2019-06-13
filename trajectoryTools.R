if(!require(lattice)) { install.packages("lattice"); library(lattice) }
if(!require(mvtnorm)) { install.packages("mvtnorm"); library(mvtnorm) }
if(!require(numDeriv)) { install.packages("numDeriv"); library(numDeriv) }
if(!require(ellipse)) { install.packages("ellipse"); library(ellipse) }
if(!require(grDevices)) { install.packages("grDevices"); library(grDevices) }
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(purrr)){install.packages("purrr"); library(purrr)}
if(!require(matrixStats)){install.packages("matrixStats"); library(matrixStats)}

mypalette = colorRampPalette(c(wes_palette(10, name = "Zissou1", type = "continuous"), "darkred"))(50)

diff_cut_off = 1e-5
max_gens = 1000
max_stand_still = 10
space_size = 10

vector_cor = function(x, y) abs(sum(Normalize(x) * Normalize(y)))

W_bar_factory = function(theta_matrix, w_cov = diag(dim(theta_matrix)[2])) {
  function(x) logSumExp(apply(theta_matrix, 1, function(theta) dmvnorm(x, mean = theta, w_cov, log = T)))
}

#start_position = rep(0, 10)
#G = G_corr
#W_bar = W_bar_multi
#W_bar(start_position)
#grad(W_bar, start_position)
#omega = diag(dim(G)[1])
#scale = 6
calculateTrajectory <- function (start_position, G, W_bar, omega = diag(dim(G)[1]), scale = 2) {
  p = dim(G)[1]
  trajectory = matrix(NA, max_gens, p)
  betas = matrix(NA, max_gens, p)
  current_position = start_position
  stand_still_counter = 0
  net_beta = rep(0, p)
  gen = 1
  while(gen <= max_gens){
    trajectory[gen,] = current_position
    beta = grad(W_bar, t(current_position))
    betas[gen,] = beta
    net_beta = net_beta + beta
    next_position = current_position + (G/scale)%*%beta
    if(Norm(next_position - current_position) < diff_cut_off){
      stand_still_counter = stand_still_counter + 1
    }
    if(stand_still_counter > max_stand_still){
      break
    }
    current_position = next_position
    gen = gen+1
  }
  trajectory = unique(trajectory[!is.na(trajectory[,1]),])
  betas = betas[!is.na(betas[,1]),]
  net_dz = trajectory[dim(trajectory)[1],] - start_position
  return(list(start_position = start_position,
              trajectory = trajectory, 
              betas = betas, 
              net_beta = net_beta,
              net_dz = net_dz))
}

plotDzgmax_normdz = function(results, gmax){
  dz_gmax = sapply(results, function(x) vector_cor(x$net_dz, gmax))
  norm_dz = sapply(results, function(x) Norm(x$net_dz))
  plot(dz_gmax, norm_dz, pch = 19, 
       xlab = expression(paste("Vector correlation between ", Delta,"z and ",g[max])),
       ylab = expression(paste("||", Delta,"z||")))
}

plotW_bar = function(W_bar, step = 0.2){
  x <- seq(-space_size, space_size, step) ## valores para mu
  y <- seq(-space_size, space_size, step)
  X <- as.matrix(expand.grid(x, y))
  Z <- vector()
  for(i in 1:nrow(X)){
    Z[i] <- W_bar(c(X[i,1], X[i,2]))
  }
  Z = exp(Z - log(sum(exp(Z))))
  b <- matrix(Z, length(x))
  mypalette = colorRampPalette(c("white", wes_palette(10, name = "Zissou1", type = "continuous"), "darkred"))
  filled.contour(x, y, z = b, color.palette = mypalette, xlim = c(-space_size, space_size), ylim = c(-space_size, space_size),
                 plot.axes = {
                   axis(1, at = seq(-10, 10));
                   axis(2, at = seq(-10, 10));
                 })
}