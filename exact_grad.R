library(mvtnorm)
library(numDeriv)
library(evolqg)
library(matrixStats)


sigma = RandomMatrix(10)
x = rep(0.9, 10)
m_1 = rep(0, 10)
m_2 = rep(1, 10)

W = function(x) logSumExp(c(dmvnorm(x, m_1, sigma, log = TRUE),
                            dmvnorm(x, m_2, sigma, log = TRUE)))
W_1 = function(x) dmvnorm(x, m_1, sigma)
W_2 = function(x) dmvnorm(x, m_2, sigma)
W_1(x)
W_2(x)
exp(W(x))
grad(W, x)

(-W_1(x) * solve(sigma, (x - m_1)) + -W_2(x) * solve(sigma, (x - m_2)))/exp(W(x))


W_bar_factory = function(theta_matrix, w_cov = diag(dim(theta_matrix)[2])) {
  function(x) logSumExp(apply(theta_matrix, 1, function(theta) dmvnorm(x, mean = theta, w_cov, log = T)))
}

W_bar_gradient_factory = function(theta_matrix, w_cov = diag(dim(theta_matrix)[2])){
 function(x) rowSums(apply(theta_matrix, 1, function(theta) - dmvnorm(x, mean = theta, w_cov) * solve(w_cov, x - theta)))/exp(W_bar_factory(theta_matrix, w_cov)(x))
}

w_cov = diag(2)
x = c(0.4, 0.5)
theta_matrix = matrix(c(0, 0, 1, 1), 2, 2, byrow = T)
W_bar = W_bar_factory(theta)
grad(W_bar, x)
W_grad = W_bar_gradient_factory(theta)
W_grad(x)
