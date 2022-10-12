options(stringsAsFactors = F)

library(umap)
library(scales)
library(MASS)
library(RColorBrewer)
library(colorRamps)
library(vegan)
library(igraph)
library(entropy)
library(jpeg)
library(repr)
library(swaRm)
library(dunn.test)
library(vioplot)
library(lsa)
library(shape)
library(gtools)
library(entropy)
library(pracma)
library(fractaldim)

###########################
##### Utility functions#####
###########################
##' split_with_overlap': Splits a vector into n chunks with x amount of overlap
split_with_overlap <- function(vec, seg_length, overlap) {
  starts <- seq(1, nrow(vec), by = seg_length - overlap)
  ends <- starts + seg_length - 1
  ends[ends > nrow(vec)] <- nrow(vec)

  lapply(1:length(starts), function(i) {
    vec[starts[i]:ends[i], ]
  })
}

##' darken_color': Darkens a given color by a defined factor amount, useful for plotting
darken_color <- function(color, factor = 1.4) {
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  col
}

##' fit_ellipse': Fits an ellipse to a given set of points in xy space
# Function taken from: http://r.789695.n4.nabble.com/Fitting-a-half-ellipse-curve-tp2719037p2720560.html
fit_ellipse <- function(x, y = NULL) {
  # Least squares fitting of an ellipse to point data using the algorithm
  # described in: Radim Halir & Jan Flusser. 1998.
  # Adapted from the original Matlab code by Michael Bedward (2010)
  # michael.bedward@gmail.com
  #
  # Subsequently improved by John Minter (2012)
  # Arguments:
  # x, y - x and y coordinates of the data points.
  #        If a single arg is provided it is assumed to be a
  #        two column matrix.
  #
  # Returns a list with the following elements:
  #
  # coef - coefficients of the ellipse as described by the general
  #        quadratic:  ax^2 + bxy + cy^2 + dx + ey + f = 0
  #
  # center - center x and y
  #
  # major - major semi-axis length
  #
  # minor - minor semi-axis length

  EPS <- 1.0e-8
  dat <- xy.coords(x, y)

  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y)
  D2 <- cbind(dat$x, dat$y, 1)
  S1 <- t(D1) %*% D1
  S2 <- t(D1) %*% D2
  S3 <- t(D2) %*% D2
  T <- -solve(S3) %*% t(S2)
  M <- S1 + S2 %*% T
  M <- rbind(M[3, ] / 2, -M[2, ], M[1, ] / 2)
  evec <- eigen(M)$vec
  cond <- 4 * evec[1, ] * evec[3, ] - evec[2, ]^2
  a1 <- evec[, which(cond > 0)]
  f <- c(a1, T %*% a1)
  names(f) <- letters[1:6]

  # calculate the center and lengths of the semi-axes
  A <-
    matrix(
      c(2 * f[1], f[2], f[2], 2 * f[3]),
      nrow = 2,
      ncol = 2,
      byrow = T
    )
  b <- matrix(c(-f[4], -f[5]),
    nrow = 2,
    ncol = 1,
    byrow = T
  )
  soln <- solve(A) %*% b
  b2 <- f[2]^2 / 4

  center <- c(soln[1], soln[2])
  names(center) <- c("x", "y")

  num <-
    2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2] * f[4] *
      f[5] / 4 - f[1] * f[3] * f[6])
  den1 <- (b2 - f[1] * f[3])
  den2 <- sqrt((f[1] - f[3])^2 + 4 * b2)
  den3 <- f[1] + f[3]

  semi_axes <-
    sqrt(c(num / (den1 * (den2 - den3)), num / (den1 * (-den2 - den3))))

  # calculate the angle of rotation
  term <- (f[1] - f[3]) / f[2]
  angle <- atan(1 / term) / 2

  list(
    coef = f,
    center = center,
    major = max(semi_axes),
    minor = min(semi_axes),
    angle = unname(angle)
  )
}

##' euc_dist': Calculates the Euclidean distance between two sets of xy
# coordinates
euc_dist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
}
