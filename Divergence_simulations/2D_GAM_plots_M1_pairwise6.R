# rm(list = ls())
library(mgcv)
library(colorRamps)
# load("300.sim.gam.results.rdata")
scripts<-"sim_scripts"
source(paste0(scripts,"/gam.plotting.funcs.R"))
op <- par(mar = c(0, 0, 0, 0))

my.pal <- function(...) {
  library(colorRamps)
  matlab.like(...)
}

n.grid <- 50

ne.lab <- expression("log"[10]*"N"[e])
div.lab <- expression("log"[10]*"T")
mut.lab <- expression("log"[10]*mu)

diag.zlim <- c(0,1)
pdf(paste0("PD.heatmaps.",".pdf"), width = 7, height = 9)

plot.mat <- matrix(1:6, ncol = 2, byrow = T)
plot.mat <- rbind(c(7, 8), plot.mat)
plot.mat <- cbind(plot.mat, c(9, 10, 10, 10))

layout(plot.mat, widths = c(3, 3, 1.5), heights = c(1.5, 10, 10, 10))

par(mar = c(5, 5.5, 1, 1) + 0.1, oma = c(0, 1, 0, 1))

gam.persp.plot(mdl1.mig.eq0.diag, "l.div.time", "l.Ne",
               xlab = div.lab, ylab = ne.lab, pal = my.pal,
               n.grid = n.grid, zlim = diag.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "a", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.gt0.diag, "l.div.time", "l.Ne",
               xlab = div.lab, ylab = ne.lab, pal = my.pal,
               n.grid = n.grid, zlim = diag.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "d", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.eq0.diag, "l.div.time", "l.mut.rate",
               xlab = div.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = diag.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "b", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.gt0.diag, "l.div.time", "l.mut.rate",
               xlab = div.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = diag.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "e", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.eq0.diag, "l.Ne", "l.mut.rate",
               xlab = ne.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = diag.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "c", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.gt0.diag, "l.Ne", "l.mut.rate",
               xlab = ne.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = diag.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "f", cex = 2, font = 2, adj = c(-0.5, 1.5))

par(mar = c(0, 5.3, 0, 1))

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, "m = 0", adj = 0.5, cex = 3, xpd = T)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, "m > 0", adj = 0.5, cex = 3, xpd = T)

plot.new()

par(mar = c(5, 3, 1, 3) + 0.1)
contour.legend(seq(0, 100, 5), my.pal(20))
axis(4, at = seq(0, 100, 10), las = 1, cex.axis = 2)
# mtext(expression("log"[10]*"d"[A]), 1, line = 2, cex = 1.5)
mtext("PD", 1, line = 2, cex = 1.5)

dev.off()
par(op)

pdf(paste0("dA.heatmaps.",".pdf"), width = 7, height = 9)

plot.mat <- matrix(1:6, ncol = 2, byrow = T)
plot.mat <- rbind(c(7, 8), plot.mat)
plot.mat <- cbind(plot.mat, c(9, 10, 10, 10))

layout(plot.mat, widths = c(3, 3, 1.5), heights = c(1.5, 10, 10, 10))

par(mar = c(5, 5.5, 1, 1) + 0.1, oma = c(0, 1, 0, 1))

gam.persp.plot(mdl1.mig.eq0.dA, "l.div.time", "l.Ne",
               xlab = div.lab, ylab = ne.lab, pal = my.pal,
               n.grid = n.grid, zlim = dA.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "a", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.gt0.dA, "l.div.time", "l.Ne",
               xlab = div.lab, ylab = ne.lab, pal = my.pal,
               n.grid = n.grid, zlim = dA.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "d", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.eq0.dA, "l.div.time", "l.mut.rate",
               xlab = div.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = dA.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "b", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.gt0.dA, "l.div.time", "l.mut.rate",
               xlab = div.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = dA.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "e", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.eq0.dA, "l.Ne", "l.mut.rate",
               xlab = ne.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = dA.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "c", cex = 2, font = 2, adj = c(-0.5, 1.5))

gam.persp.plot(mdl1.mig.gt0.dA, "l.Ne", "l.mut.rate",
               xlab = ne.lab, ylab = mut.lab, pal = my.pal,
               n.grid = n.grid, zlim = dA.zlim, line = 3.2, cex = 1.5)
text(par("usr")[1], par("usr")[4], "f", cex = 2, font = 2, adj = c(-0.5, 1.5))

par(mar = c(0, 5.3, 0, 1))

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, "m = 0", adj = 0.5, cex = 3, xpd = T)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, "m > 0", adj = 0.5, cex = 3, xpd = T)

plot.new()

par(mar = c(5, 3, 1, 3) + 0.1)
contour.legend(seq(0, 100, 5), my.pal(20))
axis(4, at = seq(0, 100, 10), labels = seq(dA.zlim[1], dA.zlim[2], ((dA.zlim[2]-dA.zlim[1])/10)), las = 1, cex.axis = 2)
mtext(expression("log"[10]*"d"[A]), 1, line = 2, cex = 1.5)

dev.off()
par(op)
