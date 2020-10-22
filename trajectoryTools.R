if(!require(lattice)) { install.packages("lattice"); library(lattice) }
if(!require(mvtnorm)) { install.packages("mvtnorm"); library(mvtnorm) }
if(!require(numDeriv)) { install.packages("numDeriv"); library(numDeriv) }
if(!require(ellipse)) { install.packages("ellipse"); library(ellipse) }
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(purrr)){install.packages("purrr"); library(purrr)}
if(!require(matrixStats)){install.packages("matrixStats"); library(matrixStats)}
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}
if(!require(NCmisc)){install.packages("NCmisc"); library(NCmisc)}
if(!require(wesanderson)){install.packages("wesanderson"); library(wesanderson)}
if(!require(cowsay)){install.packages("cowsay"); library(cowsay)}
if(!require(MASS)){install.packages("MASS"); library(MASS)}

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

randomPeaks = function(n = n_peaks, p = n_traits, x = rep(1, p), intervals = 1, prop = 1, dz_limits, 
                       max_uniform = n * 100, sigma_init = 2, sigma_step = 0.01, verbose = FALSE){
  steps = length(intervals)
  counter = vector("numeric", steps)
  n_per = ceiling(n * prop)
  peaks = matrix(0, n, p)
  k = 1
  attempts = 1
  while(k <= n & attempts < max_uniform){
    attempts = attempts + 1
    rpeak = Normalize(rnorm(p))
    corr = vector_cor(x, rpeak)
    for(i in 1:steps) {
      if(corr < intervals[i]){
        if(counter[i] < n_per[i]){
          counter[i] = counter[i] + 1
          if(verbose) print(counter)
          peaks[k,] = rpeak * runif(1, dz_limits[1], dz_limits[2])
          k = k + 1
        }
        break
      }
    }
  }
  if(k < n){
    mask = which(counter != n_per)
    mask = c(mask[1]-1, mask)
    target_intervals = intervals[mask]
    sigma = sigma_init
    while(k <= n){
      rpeak = Normalize(x + rnorm(p, 0, sigma))
      corr = vector_cor(x, rpeak)
      if(corr < target_intervals[1]){
        sigma = sigma - sigma_step
      } else if(corr > target_intervals[length(target_intervals)]) {
        sigma = sigma + sigma_step
      } else { for(i in 1:steps) {
        if(corr < intervals[i]){
          if(counter[i] < n_per[i]){
            counter[i] = counter[i] + 1
            if(verbose) print(counter)
            peaks[k,] = rpeak * runif(1, dz_limits[1], dz_limits[2])
            k = k + 1
            mask = which(counter != n_per)
            mask = c(mask[1]-1, mask)
            target_intervals = intervals[mask]
          }
          break
        }
      }
      }
    }
  }
  peaks[sample(1:n, n),]
}

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
    if(Norm(next_position) > space_size*2) stop("Out of bounds")
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
  nGen = nrow(trajectory)
  if(trim_trajectory){
      trim = round(seq(1, nGen, length.out = 100))
      trajectory = trajectory[trim,]
      betas      = betas[trim,]
  }
  return(list(start_position = start_position,
              trajectory = trajectory, 
              betas = betas, 
              net_beta = net_beta,
              net_dz = net_dz, 
              nGen = nGen))
}

runSimulation = function(G_type = c("Diagonal", "Integrated"), G = NULL,
                         n_peaks = 1, p, rho = 0.7, scale = 6, peakPool = NULL, theta = NULL){
    G_type = match.arg(G_type)
    if(is.null(G)){
        if(G_type == "Diagonal"){
          G = G_factory(p, 0.1)
          gmax = eigen(G)$vectors[,1]
        } else if(G_type == "Integrated"){
            G = G_factory(p, rho)
            gmax = eigen(G)$vectors[,1]
        } else stop("Unknown G type")
    } else
        gmax = eigen(G)$vectors[,1]
    if(n_peaks == 1){
        Surface_type = "Single"
    } else
        Surface_type = "Muliple"
    if(is.null(theta)) theta = matrix(peakPool[sample(1:nrow(peakPool), n_peaks),], n_peaks, p)
    W_bar = W_bar_factory(theta)
    W_bar_grad = W_bar_gradient_factory(theta)
    trajectory = calculateTrajectory(rep(0, p), G, W_bar, W_bar_grad, scale = scale)
    trajectory$G_type = G_type
    trajectory$G = G
    trajectory$gmax = gmax
    trajectory$theta = theta
    trajectory$z = trajectory$trajectory[dim(trajectory$trajectory)[1],]
    #trajectory$W_bar = W_bar
    #trajectory$W_bar_grad = W_bar_grad
    trajectory$Surface_type = Surface_type
    trajectory$normZ = Norm(trajectory$z)
    return(trajectory)
}

runTrypitch = function(G_diag, peakPool_diag, G_corr, peakPool_corr, n = 1000, n_peaks = 50, scale = 4, parallel = TRUE){
  say("G diag single")
  G_diag_W_single = llply(1:n, function(x) runSimulation("Diag", G_diag,   1, n_traits, 
                                                         scale = scale, 
                                                         peakPool = peakPool_diag), 
                          .parallel = parallel)
  say("G corr single")
  G_corr_W_single = llply(1:n, function(x) runSimulation("Inte", G_corr,   1, n_traits, 
                                                         scale = scale, 
                                                         peakPool = peakPool_corr), 
                          .parallel = parallel)
  say("G diag multi")
  G_diag_W_multi  = llply(1:n, function(x) runSimulation("Diag", G_diag, n_peaks, n_traits, 
                                                         scale = scale, 
                                                         peakPool = peakPool_diag), 
                          .parallel = parallel)
  say("G corr multi")
  G_corr_W_multi  = llply(1:n, function(x) runSimulation("Inte", G_corr, n_peaks, n_traits, 
                                                         scale = scale, 
                                                         peakPool = peakPool_corr), 
                          .parallel = parallel)
  list(DS = G_diag_W_single,
       CS = G_corr_W_single,                                                                  
       DM = G_diag_W_multi,
       CM = G_corr_W_multi)
}

runTrypitchList = function(Gs, peakPools, n_peaks, n = 1000, scale = 4, parallel = TRUE){
    n_sim = 2*length(Gs)
    results = vector("list", n_sim)
    for(i in seq(1, n_sim, 2))
    {
        results[[i]] = llply(1:n, function(x) runSimulation("Integrated", Gs[[i]],   1, n_traits, 
                                                            scale = scale, 
                                                            peakPool = peakPools[[i]]), 
                             .parallel = parallel)  
        results[[i+1]] = llply(1:n, function(x) runSimulation("Integrated", Gs[[i]], n_peaks, n_traits, 
                                                              scale = scale, 
                                                              peakPool = peakPools[[i]]), 
                               .parallel = parallel)  
    }
    names(results) = names(Gs)
    return(results)
}

ReplaceDiagonal = function(x, d){
    d = sqrt(d)
    c.x = cov2cor(x)
    outer(d, d) * c.x
}
make_matrix = function(eVal, eVec, p = 1) eVec %*% diag(eVal^p) %*% t(eVec)
expEigenVal = function(mat, p){
  eigX = eigen(mat)
  eVal = eigX$values
  eVec = eigX$vectors
  new_mat = ReplaceDiagonal(make_matrix(eVal, eVec, p), d = diag(mat))
  return(new_mat)
}

G_factory = function(p, rho, sigma = 0.1){
  while(TRUE){
    G = matrix(rnorm(p*p, rho, sigma), p, p)
    G = (G + t(G))/2
    diag(G) = rnorm(p, 1, sigma)
    tryCatch({chol(G); break}, error = function(x) FALSE)
  }
  G
}
