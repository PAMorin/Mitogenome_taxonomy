
sim.gam <- function(form, df, label) {
  library(colorRamps)
  df$right <- round(df$sample.size * 2 * (df$diagnosability/100), 0)
  df$wrong <- df$sample.size * 2 - df$right
  gam(form, data = df, family = binomial)
}


gam.xyz <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, type = "link", ...) {
  library(mgcv)
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      } else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) break
    }
    if (k < 2) stop("Model does not seem to have enough terms to do anything useful")
  } else {
    if (sum(view %in% v.names) != 2) 
      stop(paste(c("view variables must be one of", v.names), collapse = ", "))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  } else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(paste("View variables must contain more than one value. view = c(", view[1], ",", view[2], ").", sep = ""))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  
  z <- matrix(fv$fit, n.grid, n.grid)    
  list(m1 = m1, m2 = m2, z = z)
}


my.filled.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, 
                               xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                               levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                               col = color.palette(length(levels) - 1)) 
{
  plot.new()
  plot.window(xlim, ylim, "", xaxs = "i", yaxs = "i", asp = NA)
  .filled.contour(x, y, z, levels, col)
}


gam.persp.plot <- function(gam.fit, x.var, y.var, xlab = NULL, ylab = NULL, 
                           pal = rainbow, n.grid = 100, zlim = c(0,1), line, cex, ...) 
{
  xyz <- gam.xyz(gam.fit, view = c(x.var, y.var), n.grid = n.grid, type = "response")
  with(xyz, my.filled.contour(m1, m2, z, zlim = zlim, color.palette = pal, ...))
  axis(1, lwd = 2, cex.axis = cex)
  axis(2, lwd = 2, cex.axis = cex)
  box(lwd = 2)
  mtext(xlab, 1, line = line, cex = cex)
  mtext(ylab, 2, line = line, cex = cex)
}

contour.legend <- function(levels, col, cex = 1) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  box()
}
