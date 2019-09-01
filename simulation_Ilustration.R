source("./trajectoryTools.R")

space_size = 10
npeaks = 100
p = 2

uniform_peaks = matrix(runif(p*npeaks, -10, 10), npeaks, p)
plot(uniform_peaks, axes = FALSE, pch = 17, xlab = "", ylab= "", xlim=c(-10.1, 10.1), ylim = c(-10.1, 10.1), asp = 1)
polygon(x = matrix(c(-10, -10, 
                     -10, 10,
                     10, 10, 
                     10, -10,
                     -10, -10), 5, 2, byrow = T))
abline(v=0)
abline(h=0)

peaks = randomPeaks(npeaks, p, intervals = c(1), prop = c(1), dz_lim = c(0, 10))
plot(peaks, axes = FALSE, pch = 17, xlab = "", ylab= "", xlim=c(-10.1, 10.1), ylim = c(-10.1, 10.1), asp = 1)
draw.circle(0, 0, 10, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
abline(v=0)
abline(h=0)

peaks = randomPeaks(npeaks, p, intervals = c(1), prop = c(1), dz_lim = c(3, 10))
plot(peaks, axes = FALSE, pch = 17, xlab = "", ylab= "", xlim=c(-10.1, 10.1), ylim = c(-10.1, 10.1), asp = 1)
draw.circle(0, 0, 10, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
draw.circle(0, 0, 3, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
abline(v=0)
abline(h=0)
gplotW_bar(peaks, 10)

space_size = 10
npeaks = 50
p = 2
set.seed(50) # 50
{
  peakPool = randomPeaks(100, p = p, dz_limits = c(3, 6), 
                         intervals = c(1), prop = c(1))
  theta = matrix(peakPool[sample(1:nrow(peakPool), npeaks),], npeaks, p)
  x = runSimulation("Integrated", rho = 0.75, p = p, scale = 4, theta = theta)
  p_x = gplotW_bar_trajectory(x, 9, main = "Integrated - rho = 0.75") 
  y = runSimulation("Diagonal", p = 2, scale = 4, theta = theta)
  p_y = gplotW_bar_trajectory(y, 9, main = "Diagonal") 
  pxy = plot_grid(p_x , p_y)

}
save_plot("~/surfaces.png", pxy, base_height = 7, base_asp = 1, ncol = 2)

space_size = 6
npeaks = 5
set.seed(2) # 2, 10, 15, 18, 22, 25, 27, 30, 31, 38, 45, 47, 48, 49
{
  peakPool = randomPeaks(100, p = 2, dz_limits = c(3, space_size), 
                         intervals = c(1), prop = c(1))
  theta = matrix(peakPool[sample(1:nrow(peakPool), npeaks),], npeaks, p)
  #plotW_bar(theta, mypalette = colorRampPalette(wes_palette(name = "Zissou1", type = "continuous")))
  x = runSimulation("Integrated", rho = 0.9, n_peaks = 5, p = 2, scale = 4, theta = theta)
  p_x = gplotW_bar_trajectory(x, 8)
  y = runSimulation("Diagonal", n_peaks = 5, p = 2, scale = 4, theta = theta)
  p_y = gplotW_bar_trajectory(y, 8)
  plot_grid(p_x, p_y, labels = c("Integrated", "Diagonal"))
}
