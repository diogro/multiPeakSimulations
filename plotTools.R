if(!require(remotes)){install.packages("remotes"); library(remotes)}
if(!require(ContourFunctions)){remotes::install_github("CollinErickson/contour"); library(ContourFunctions)}
if(!require(plotrix)){install.packages("plotrix"); library(plotrix)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
theme_set(theme_cowplot())
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }
if(!require(grDevices)) { install.packages("grDevices"); library(grDevices) }

plotDzgmax_normdz = function(results, ylim, main = ""){
  df = data.frame(dz_gmax = sapply(results, function(x) vector_cor(x$net_dz, x$gmax)),
                  norm_dz = sapply(results, function(x) Norm(x$net_dz)))
  ggplot(df, aes(dz_gmax, norm_dz)) + geom_point(shape = 19) +  
    scale_y_continuous(limits = ylim) + scale_x_continuous(limits = c(0, 1)) +
    labs(x = expression(paste("Vector correlation between ", Delta,"z and ",g[max])),
         y = expression(paste("||", Delta,"z||"))) + 
    ggtitle(main)
}

plotW_bar_trajectory = function(run, xlimits = c(-space_size, space_size), ylimits = c(-space_size, space_size),
                                mypalette = colorRampPalette(c("white", wes_palette(10, name = "Zissou1", type = "continuous"), "darkred")), 
                                main = "", ...){
  W_bar = W_bar_factory(run$theta)
  x <- seq(xlimits[1], xlimits[2], resolution) 
  y <- seq(ylimits[1], ylimits[2], resolution) 
  Z <- matrix(NA, length(x), length(y))
  for(i in 1:length(x))
    for(j in 1:length(y))
      Z[i, j] <- W_bar(c(x[i], y[j]))
  if(log) { Z = Z - logSumExp(Z)
  } else Z = exp(Z - logSumExp(Z))
  cf_grid(x = x, y = y, z = Z, xlim = xlimits, ylim = ylimits, color.palette = mypalette, main = main,
          afterplotfunc=function() {
            points(run$theta, pch = 17)
            points(run$trajectory, pch = 19)
            abline(v=0)
            abline(h=0)
          }, ...)
}

plotW_bar = function(theta, space_size = 6, xlimits = c(-space_size, space_size), ylimits = c(-space_size, space_size), resolution = 0.2, 
                     mypalette = colorRampPalette(c("white", wes_palette(10, name = "Zissou1", type = "continuous"), "darkred")), 
                     main = "", log = FALSE, ...){
  W_bar = W_bar_factory(theta)
  x <- seq(xlimits[1], xlimits[2], resolution) 
  y <- seq(ylimits[1], ylimits[2], resolution) 
  Z <- matrix(NA, length(x), length(y))
  for(i in 1:length(x))
    for(j in 1:length(y))
      Z[i, j] <- W_bar(c(x[i], y[j]))
  if(log) { Z = Z - logSumExp(Z)
  } else Z = exp(Z - logSumExp(Z))
  cf_grid(x = x, y = y, z = Z, xlim = xlimits, ylim = ylimits, color.palette = mypalette, main = main, ...,
          afterplotfunc=function() {
            points(theta, pch = 17)
            abline(v=0)
            abline(h=0)
          })
}

gplotW_bar = function(theta, space_size = 6, xlimits = c(-space_size, space_size), ylimits = c(-space_size, space_size), resolution = 0.2,
                      mypalette = colorRampPalette(c("white", wes_palette(10, name = "Zissou1", type = "continuous"), "darkred")), 
                      log = FALSE, main = "", ...){
  W_bar = W_bar_factory(theta)
  x <- seq(xlimits[1], xlimits[2], resolution) 
  y <- seq(ylimits[1], ylimits[2], resolution) 
  X <- as.matrix(expand.grid(x, y))
  Z <- vector()
  for(i in 1:nrow(X)){
    Z[i] <- W_bar(c(X[i,1], X[i,2]))
  }
  if(log) { Z = Z - logSumExp(Z)
  } else Z = exp(Z - logSumExp(Z))
  gcf_grid(x, y, Z, xlim = xlimits, ylim = ylimits, color.palette = mypalette, 
           main = main, mainminmax = FALSE, mainminmax_minmax = FALSE, ...) + 
    geom_point(data=data.frame(theta), aes(X1, X2), shape = 17) + ggtitle(main) + 
    coord_fixed() + theme_void() + theme(legend.position = "none") +
    geom_segment(aes(x = 0, xend = 0, y = ylimits[1], yend = ylimits[2])) + geom_segment(aes(y = 0, yend = 0, x = xlimits[1], xend = xlimits[2]))
}

gplotW_bar_trajectory = function(run, space_size = 6, xlimits = c(-space_size, space_size), ylimits = c(-space_size, space_size), resolution = 0.2,
                                 mypalette = colorRampPalette(c("white", wes_palette(10, name = "Zissou1", type = "continuous"), "darkred")), 
                                 log = FALSE, main = "", ...){
  W_bar = W_bar_factory(run$theta)
  x <- seq(xlimits[1], xlimits[2], resolution) 
  y <- seq(ylimits[1], ylimits[2], resolution) 
  X <- as.matrix(expand.grid(x, y))
  Z <- vector()
  for(i in 1:nrow(X)){
    Z[i] <- W_bar(c(X[i,1], X[i,2]))
  }
  if(log) { Z = Z - logSumExp(Z)
  } else Z = exp(Z - logSumExp(Z))
  gcf_grid(x, y, Z, xlim = xlimits, ylim = ylimits, color.palette = mypalette, mainminmax = FALSE, mainminmax_minmax = FALSE, ...) + 
    geom_point(data=data.frame(run$theta), aes(X1, X2), shape = 17) + 
    geom_point(data=data.frame(run$trajectory), aes(X1, X2), shape = 19) + 
    ggtitle(main) + coord_fixed() + theme_void() + theme(legend.position = "none") +
    geom_segment(aes(x = 0, xend = 0, y = ylimits[1], yend = ylimits[2])) + geom_segment(aes(y = 0, yend = 0, x = xlimits[1], xend = xlimits[2]))
}
