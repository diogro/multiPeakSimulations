source("./trajectoryTools.R")

npeaks = 50
p = 2

uniform_peaks = matrix(runif(p*npeaks, -6, 6), npeaks, p)
plot(uniform_peaks, axes = FALSE, pch = 17, xlab = "", ylab= "", xlim=c(-6.1, 6.1), ylim = c(-6.1, 6.1), asp = 1)
polygon(x = matrix(c(-6, -6, 
                     -6, 6,
                     6, 6, 
                     6, -6,
                     -6, -6), 5, 2, byrow = T))
abline(v=0)
abline(h=0)

peaks = randomPeaks(npeaks, p, intervals = c(1), prop = c(1), dz_lim = c(0, 6))
plot(peaks, axes = FALSE, pch = 17, xlab = "", ylab= "", xlim=c(-6.1, 6.1), ylim = c(-6.1, 6.1), asp = 1)
draw.circle(0, 0, 6, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
abline(v=0)
abline(h=0)

peaks = randomPeaks(npeaks, p, intervals = c(1), prop = c(1), dz_lim = c(3, 5.9))
plot(peaks, axes = FALSE, pch = 17, xlab = "", ylab= "", xlim=c(-6.1, 6.1), ylim = c(-6.1, 6.1), asp = 1)
draw.circle(0, 0, 6, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
draw.circle(0, 0, 3, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
abline(v=0)
abline(h=0)



peakPool = 
plot(apply(peakPool[sample(1:10000, 200),], 1, function(x) vector_cor(x, rep(1, 8))))
hist((apply(peakPool_G_diag, 1, Norm)))
prcomp(peakPool)

peaks = randomPeaks(5000, 8, intervals = c(0.7, 0.8, 1), prop = c(0.9, 0.06, 0.04), dz_lim = c(3, 6))
hist(sort(apply(peaks, 1, vector_cor, rep(1, 8))), breaks = 50)

set.seed(32) # 2, 10, 15, 18, 22, 25, 27, 30, 31
{
par(mfrow=c(1,1))
peakPool = randomPeaks(100, p = 2, dz_limits = c(3, 6), 
                                         intervals = c(1), prop = c(1))
theta = matrix(peakPool[sample(1:nrow(peakPool), 5),], 5, p)
#plotW_bar(theta, mypalette = colorRampPalette(wes_palette(name = "Zissou1", type = "continuous")))
x = runSimulation("Integrated", rho = 0.9, n_peaks = 5, p = 2, scale = 4, theta = theta)
plotW_bar_trajectory(x, asp = 1)
y = runSimulation("Diagonal", n_peaks = 5, p = 2, scale = 4, theta = theta)
plotW_bar_trajectory(y, asp =1)
}
