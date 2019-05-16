if(!require(lattice)) { install.packages("lattice"); library(lattice) }
if(!require(mvtnorm)) { install.packages("mvtnorm"); library(mvtnorm) }
if(!require(numDeriv)) { install.packages("numDeriv"); library(numDeriv) }
if(!require(ellipse)) { install.packages("ellipse"); library(ellipse) }
if(!require(grDevices)) { install.packages("grDevices"); library(grDevices) }
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}

diff_cut_off = 1e-4
max_gens = 100
max_stand_still = 10

calculateTrajectory <- function (G, W_bar, omega = diag(dim(G)[1]), start_position = rep(0, dim(G)[1]), scale = 2) {
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
      trajectory = unique(trajectory[!is.na(trajectory[,1]),])
      betas = betas[!is.na(betas[,1]),]
      break
    }
    net_dz = trajectory[dim(trajectory)[1],] - start_position
    current_position = next_position
    gen = gen+1
  }
  return(list(trajectory = trajectory, 
              betas = betas, 
              net_beta = net_beta,
              net_dz = net_dz))
}

w_cov = matrix(c(1.0, 0.0,
                 0.0, 1.0), ncol = 2)


W_bar_multi = function(x) {
  log(
    dmvnorm(x, mean = c(4, 5), sigma = w_cov) +
      1.1*dmvnorm(x, mean = c(1, 5), sigma = w_cov) +
      dmvnorm(x, mean = c(2, 2), sigma = w_cov) +
      dmvnorm(x, mean = c(7, 5), sigma = w_cov) +
      dmvnorm(x, mean = c(5, 2), sigma = w_cov))
}

W_bar_single = function(x) {
  log(dmvnorm(x, mean = c(3, 3), sigma = w_cov))
}

G_corr = matrix(c(1.1, 0.8,
             0.8, 1.0), ncol = 2)

G_diag = matrix(c(1.1, 0.0,
                  0.0, 1.0), ncol = 2)

n_sims = 100
results = vector("list", n_sims)
for(i in 1:n_sims){
  random_start = runif(2, -10, 10)
  results[[i]] = calculateTrajectory(G, W_bar, start_position = random_start, scale = 10)
}

mypalette = colorRampPalette(c(wes_palette(10, name = "Zissou1", type = "continuous"), "darkred"))(n_sims)

plot(results[[1]]$trajectory, xlim = c(-10, 10), ylim = c(-10, 10))
for(i in 2:n_sims){
  points(results[[i]]$trajectory, col=mypalette[i])  
}
