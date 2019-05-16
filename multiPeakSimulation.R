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

calculateTrajectory <- function (G, W, omega = diag(dim(G)[1]), start_position = rep(0, dim(G)[1]), scale = 2) {
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
out = calculateTrajectory(G, NULL, start_position = c(1.5, 0), scale = 10)
plot(out$trajectory)
