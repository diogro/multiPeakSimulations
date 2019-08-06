if(!require(lattice)) { install.packages("lattice"); library(lattice) }
if(!require(mvtnorm)) { install.packages("mvtnorm"); library(mvtnorm) }
if(!require(numDeriv)) { install.packages("numDeriv"); library(numDeriv) }
if(!require(ellipse)) { install.packages("ellipse"); library(ellipse) }
if(!require(grDevices)) { install.packages("grDevices"); library(grDevices) }
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(purrr)){install.packages("purrr"); library(purrr)}
if(!require(matrixStats)){install.packages("matrixStats"); library(matrixStats)}
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}

mypalette = colorRampPalette(c(wes_palette(10, name = "Zissou1", type = "continuous"), "darkred"))(50)

vector_cor = function(x, y) abs(x %*% y/(Norm(x)*Norm(y)))

W_bar_factory = function(theta_matrix, w_cov = diag(dim(theta_matrix)[2])) {
  function(x) logSumExp(apply(theta_matrix, 1, function(theta) dmvnorm(x, mean = theta, w_cov, log = T)))
}

W_bar_gradient_factory = function(theta_matrix, w_cov = NULL){
    if(is.null(w_cov)){
  function(x) rowSums(apply(theta_matrix, 1, function(theta) - dmvnorm(x, mean = theta) * t(x - theta)))/exp(W_bar_factory(theta_matrix)(x))
    } else{
  function(x) rowSums(apply(theta_matrix, 1, function(theta) - dmvnorm(x, mean = theta, w_cov) * solve(w_cov, x - theta)))/exp(W_bar_factory(theta_matrix, w_cov)(x))
    }
}

randomPeaks = function(n = n_peaks, p = n_traits, x = rep(1, p), steps = 10, minimum, maximum){
  counter = vector("numeric", steps)
  n_per = n/steps
  corr_step = 1/steps
  peaks = matrix(0, n, p)
  k = 1
  while(k <= n){
    rpeak = Normalize(rnorm(p))
    corr = vector_cor(x, rpeak)
    for(i in 1:steps) {
      if(corr < i*corr_step){
        if(counter[i] < n_per){
          counter[i] = counter[i] + 1
          peaks[k,] = rpeak * runif(1, minimum, maximum)
          k = k + 1
          break
        }
        break
      }
    }
  }
  peaks[sample(1:n, n),]
}

#start_position = rep(0, n_traits)
#G = G_corr
#W_bar = W_bar_single
#W_bar_single_grad = W_bar_gradient_factory(theta_single)
#W_bar_gradient = W_bar_single_grad
#W_bar(start_position)
#grad(W_bar, start_position)
#W_bar_gradient(as.vector(start_position))
#scale = 6
calculateTrajectory <- function (start_position, G, W_bar, W_bar_grad, scale = 2) {
  p = dim(G)[1]
  trajectory = matrix(NA, max_gens, p)
  betas = matrix(NA, max_gens, p)
  current_position = start_position
  stand_still_counter = 0
  net_beta = rep(0, p)
  gen = 1
  while(gen <= max_gens){
    trajectory[gen,] = current_position
    beta = W_bar_grad(as.vector(current_position))
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

plotDzgmax_normdz = function(results, ...){
  dz_gmax = sapply(results, function(x) vector_cor(x$net_dz, x$gmax))
  norm_dz = sapply(results, function(x) Norm(x$net_dz))
  plot(dz_gmax, norm_dz, pch = 19, 
       xlab = expression(paste("Vector correlation between ", Delta,"z and ",g[max])),
       ylab = expression(paste("||", Delta,"z||")), ...)
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
